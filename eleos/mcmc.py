import emcee
import schwimmbad
import copy
import sys
import os
import uuid
import time
import shutil
import numpy as np
from pathlib import Path
from mpi4py import MPI

from . import cores


class Parameter:
    def __init__(self, profile_name, param_name):
        self.profile_name = profile_name
        self.param_name = param_name

    def __str__(self):
        return f"<Parameter: {self.profile_name} {self.param_name}>"
    __repr__ = __str__


def find_variable_parameters(profiles, pct_threshold=1e-5):
    out = []
    for profile in profiles:
        params = profile.CONSTANTS + profile.VARIABLES
        for param in params:
            try:
                val = getattr(profile, param)
                err = getattr(profile, f"{param}_error")
                if err / val > pct_threshold:
                    out.append(Parameter(profile.label, param))
            except AttributeError:
                continue
    return out


def get_unique_id():
    return f"{os.getpid()}-{int(time.time()*1e9)}-{uuid.uuid4().hex[:8]}"


def log_probability(state_vector, 
                    template_core, 
                    parent_directory, 
                    params_to_fit,
                    save_below_threshold):
        
        # copy the template core
        core = copy.deepcopy(template_core)
        core.forward = True

        # set the new directory
        core.id_ = get_unique_id()
        print(f"Core ID generated: {core.id_}")
        core.parent_directory = parent_directory
        core.directory = core.parent_directory / f"core_{core.id_}"
        
        # unpack params into core
        parameters = dict(zip(params_to_fit, state_vector))
        for param, value in parameters.items():
            setattr(core.profiles[param.profile_name], param.param_name, value)
            
        # run the core with those params
        try:
            print(f"Core {core.id_} starting run")
            res = core.run(confirm=False)
        except Exception as e:
            # if the core fails then assume the params were so far off the fit failed
            print(f"Core {core.id_} failed to fit: Exception: {e}")
            core.delete_core_directory(confirm=False)
            return -np.inf
            
        # check for nan chi-sq
        if np.isnan(res.chi_sq):
            print(f"Core {core.id_} finished with a chi-sq of NaN after {res.elapsed_time}")
            core.delete_core_directory(confirm=False)
            return -np.inf

        # Print chi-sq value of the iteration
        print(f"Core {core.id_} finished with a chi-sq of {res.chi_sq} after {res.elapsed_time}")

        # Save the core if chi-squared below threshold
        if save_below_threshold is not None:
            if res.chi_sq < save_below_threshold:
                print(f"Core {core.id_} chi-sq {res.chi_sq} below threshold {save_below_threshold}, saving core")
                shutil.copytree(core.directory, parent_directory / "below_threshold" / f"core_{core.id_}")

        # Reduce the core filesize or delete the core entirely
        #core.delete_auxillary_files()
        core.delete_core_directory(confirm=False)

        # return log chisq
        return np.log(res.chi_sq)


def run_mcmc(parent_directory, 
             template_core, 
             params_to_fit, 
             nwalkers, 
             max_steps,
             resume=False,
             save_below_threshold=None):
    
    # Make directory tree if required
    parent_directory = Path(parent_directory)
    if save_below_threshold is not None:
        os.makedirs(parent_directory / "below_threshold", exist_ok=True)

    # Get the starting state vectors for each walker
    state_vectors = []
    for i in range(nwalkers):
        state_vectors.append([])
        for param in params_to_fit:
            # Get the value of the parameter and its error
            value = getattr(template_core.profiles[param.profile_name], param.param_name)
            error = getattr(template_core.profiles[param.profile_name], f"{param.param_name}_error")

            # Use the value and error as mean and std for sampling from a normal distribution
            x = np.random.normal(value, error)
            x = np.clip(x, 0, np.inf)
            state_vectors[-1].append(x)

    # set up MPI pool
    pool = schwimmbad.MPIPool(use_dill=True)

    # Set up the backend
    # Don't forget to clear it in case the file already exists
    filename = parent_directory / "sampler_chain.hdf5"
    backend = emcee.backends.HDFBackend(filename)
    if not resume:
        backend.reset(nwalkers, len(params_to_fit))

    # Set up the sampler
    sampler = emcee.EnsembleSampler(nwalkers, 
                                    len(params_to_fit), 
                                    log_probability, 
                                    kwargs={"template_core": template_core, 
                                            "parent_directory": parent_directory,
                                            "params_to_fit": params_to_fit,
                                            "save_below_threshold": save_below_threshold},
                                    backend=backend,
                                    pool=pool)   

    # Run the MCMC:
    # We'll track how the average autocorrelation time estimate changes
    index = 0
    autocorr = np.empty(max_steps)
    old_tau = np.inf

    # Now we'll sample for up to max_n steps
    for sample in sampler.sample(state_vectors, iterations=max_steps, progress=True):
        # Only check convergence every 10 steps
        if sampler.iteration % 10:
            continue

        # Compute the autocorrelation time so far
        tau = sampler.get_autocorr_time(tol=0)
        autocorr[index] = np.mean(tau)
        index += 1

        # Check convergence - change of autocorrelation time < 1% is converged
        converged = np.all(tau * 100 < sampler.iteration)
        converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
        if converged:
            break
        old_tau = tau

    # Delete the kwargs to the log_probability function
    # this allows pickling of the sampler object
    sampler.log_prob_fn.kwargs = dict()

    # close the pool
    pool.close()

    return sampler


def get_parameter_names(params_to_fit):
    names = []
    for param in params_to_fit:
        names.append(str(param))
    return names

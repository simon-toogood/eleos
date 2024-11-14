# Eleos - A Python wrapper for NEMESIS


## Current hierarchy

* constants.py
  * GASES: lookup table for gas ID/names
  * DISTANCES: planetary orbit radii in AU
* main.py
  * A scratch file for testing functionality
* profiles.py
  * Profile object: Base profile object
    * TemperatureProfile
    * GasProfile
    * Aerosol Profile
  * read_profile_string: Convert a NEMESIS string (eg 20 0 1) to a Profile object
* shapes.py
  * Shape object: Base Shape object
    * Shape0
    * Shape1
    * Shape32
    * to be expanded as needed
  * get_shape_from_id: Take a NEMESIS profile shape ID and return the corresponding Shape subclass
  * ALL_SHAPES: list of all Shape subclasses
* results.py
  * NEMESISResult: Object for storing the results of a NEMESIS core
  * \



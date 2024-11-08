import constants


def read_profile_string(string):
    tokens = [int(x) for x in string.split()]
    if tokens == [0, 0, 0]:
        return TemperatureProfile()
    elif tokens[0] > 0:
        return GasProfile(*tokens)
    else:
        raise NotImplementedError


class Profile:
    def __init__(self):
        pass

    def add_result(self, result):
        self.result = result


class TemperatureProfile(Profile):
    def __init__(self):
        super().__init__()

    def __str__(self):
        return f"TemperatureProfile [{self.create_nemesis_string()}]"

    def create_nemesis_string(self):
        return f"0 0 0"
    
    def create_apr_data(self):
        pass


class GasProfile(Profile):
    def __init__(self, gas_id, isotope_id=0, profile_shape_id=1):
        super().__init__()
        self.gas_id = gas_id
        self.isotope_id = isotope_id
        self.profile_shape_id = profile_shape_id
        self.gas_name = constants.GASES.loc[constants.GASES.radtrans_id == gas_id].name

    @classmethod
    def from_string(cls, string):
        return cls(*string.split())
    
    def __str__(self):
        return f"GasProfile({self.name}) [{self.create_nemesis_string()}]"

    def create_nemesis_string(self):
        return f"{self.gas_id} {self.isotope_id} {self.profile_shape_id}"
    
    def create_apr_data(self, *args):
        pass


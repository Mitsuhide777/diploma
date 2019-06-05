class WrongParameterName(Exception):
    def __init__(self, message):
        self.message = message


class DistributionNotSupported(Exception):
    def __init__(self, message):
        self.message = message


class EmptyDistributionsDict(Exception):
    def __init__(self, message):
        self.message = message


class EmptyParametersDict(Exception):
    def __init__(self, message):
        self.message = message


class EmptyParameterList(Exception):
    def __init__(self, message):
        self.message = message


class FileParseException(Exception):
    def __init__(self, message):
        self.message = message

class AlvPossibleFormatError(Exception):
    '''
    Cannot recognize the input format
    '''
    def __init__(self, expression):
        self.expression = 'HUBBA'
        self.message = 'File format not recognized.'

class AlvEmptyAlignment(Exception):
    def __init__(self):
        self.message = 'Input alignment contains no sequences.'

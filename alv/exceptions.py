class AlvPossibleFormatError(Exception):
    '''
    When reading a file yields no records, so probably wrong format.
    '''
    def __init__(self, expression, filename):
        self.expression = 'HUBBA'
        self.message = 'No sequence records in ' + filename

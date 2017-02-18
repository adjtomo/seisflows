
class ParameterError(ValueError):
    def __init__(self, *args):
        if len(args) == 0:
            msg = 'Bad parameter.'
            super(ParameterError, self).__init__(msg)
        elif len(args) == 1:
            msg = 'Bad parameter: %s' % args[0]
            super(ParameterError, self).__init__(msg)
        elif args[1] not in args[0]:
            msg = '%s is not defined.' % args[1]
            super(ParameterError, self).__init__(msg)
        elif key in obj:
            msg = '%s has bad value: ' % args[0], args[1].__getattr__(args[0])
            super(ParameterError, self).__init__(msg)


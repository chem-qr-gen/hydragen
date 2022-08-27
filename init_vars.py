from inspect import getargspec
def init_vars(func):
    def returnFunc(*args, **kwargs):
        selfVar = args[0]

        argSpec = getargspec(func)
        argumentNames = argSpec[0][1:]
        defaults = argSpec[3]
        if defaults is not None:
            defaultArgDict = dict(zip(reversed(argumentNames), reversed(defaults)))
            selfVar.__dict__.update(defaultArgDict)

        argDict = dict(zip(argumentNames, args[1:]))
        selfVar.__dict__.update(argDict)


        validKeywords = set(kwargs) & set(argumentNames)
        kwargDict = {k: kwargs[k] for k in validKeywords}
        selfVar.__dict__.update(kwargDict)

        func(*args, **kwargs)

    return returnFunc
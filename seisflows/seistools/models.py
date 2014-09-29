
def loadascii(dir):
    wildcard = os.path.join(dir,'*.ascii')
    model = {}
    for file in glob.glob(wildcard):
      key = os.path.splitext(os.path.basename(file))[0]
      model[key] = np.loadtxt(file)
    return model


def loadsep():
    pass


def loadnc():
    pass



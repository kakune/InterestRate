import configparser


class Parameters:
    def __init__(self):
        self.mData = {}
        self.mNameCurrentSection = None
        self.mNameCommonSection = "COMMON"

    def readParameters(self, inPathFile):
        self.mConfig = configparser.ConfigParser()
        self.mConfig.read(inPathFile)

    def __getitem__(self, inNameSection):
        return self.mConfig[inNameSection]

    def __call__(self, inNameParameter):
        lStr = None
        if (self.mConfig.has_option(self.mNameCurrentSection, inNameParameter)):
            lStr = self.mConfig[self.mNameCurrentSection][inNameParameter]
        else:
            lStr = self.mConfig[self.mNameCommonSection][inNameParameter]
        try:
            return float(lStr)
        except ValueError:
            return lStr

    def setNameCommonSection(self, inNameCommonSection):
        self.mNameCommonSection = inNameCommonSection

    def setNameCurrentSection(self, inNameCurrentSection):
        self.mNameCurrentSection = inNameCurrentSection

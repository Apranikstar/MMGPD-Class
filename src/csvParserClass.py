import csv
### HOW to run: profileFunctionParameters("HGAG23", "H", "Set9").get_flavour_values("dbar")
class getProfileFunctionParameters:
    def __init__(self,  analysisType=None, gpdType=None, analysisSet=None):
        self.dataFilename = "src/data"
        self.analysisType = analysisType
        self.gpdType = gpdType
        self.analysisSet = analysisSet

    def __get_flavour_values__(self, flavorInput):

        with open(self.dataFilename+"/"+self.analysisType+"/"+self.gpdType+"/"+self.analysisSet+".csv", newline='') as csvfile:
            csvreader = csv.reader(csvfile)
            for row in csvreader:
                flavour = row[0]
                if flavour == flavorInput:
                    return [eval(value) for value in row[1:]]

    def __call__(self, flavor):
        return self.__get_flavour_values__(flavor)


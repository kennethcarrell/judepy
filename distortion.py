import re
import statistics

####
#### a class to read in distortion information from Astrometrica ouput
####
class Distortion:
    ## define variables
    _verbose_ = False

    N_fields = 0
    x_eqn = ''
    y_eqn = ''
    N_stars = []

    x_const = 0.
    x_const_err = 0.
    x_const_arr = []
    x_xp = 0.
    x_xp_err = 0.
    x_xp_arr = []
    x_yp = 0.
    x_yp_err = 0.
    x_yp_arr = []
    x_xp2 = 0.
    x_xp2_err = 0.
    x_xp2_arr = []
    x_xpyp = 0.
    x_xpyp_err = 0.
    x_xpyp_arr = []
    x_yp2 = 0.
    x_yp2_err = 0.
    x_yp2_arr = []
    x_xp3 = 0.
    x_xp3_err = 0.
    x_xp3_arr = []
    x_xp2yp = 0.
    x_xp2yp_err = 0.
    x_xp2yp_arr = []
    x_xpyp2 = 0.
    x_xpyp2_err = 0.
    x_xpyp2_arr = []
    x_yp3 = 0.
    x_yp3_err = 0.
    x_yp3_arr = []

    y_const = 0.
    y_const_err = 0.
    y_const_arr = []
    y_xp = 0.
    y_xp_err = 0.
    y_xp_arr = []
    y_yp = 0.
    y_yp_err = 0.
    y_yp_arr = []
    y_xp2 = 0.
    y_xp2_err = 0.
    y_xp2_arr = []
    y_xpyp = 0.
    y_xpyp_err = 0.
    y_xpyp_arr = []
    y_yp2 = 0.
    y_yp2_err = 0.
    y_yp2_arr = []
    y_xp3 = 0.
    y_xp3_err = 0.
    y_xp3_arr = []
    y_xp2yp = 0.
    y_xp2yp_err = 0.
    y_xp2yp_arr = []
    y_xpyp2 = 0.
    y_xpyp2_err = 0.
    y_xpyp2_arr = []
    y_yp3 = 0.
    y_yp3_err = 0.
    y_yp3_arr = []

    ## class constructor
    def __init__(self,fname='',vb=False,reset=True):
        self._verbose_ = vb
        if(fname != ''):
            self.readAstrometrica(fname,reset)

    ## give variables some initial value
    def clearVars(self):
        self.N_fields = 0
        self.x_eqn = ''
        self.y_eqn = ''
        self.N_stars = []

        self.x_const = 0.
        self.x_const_err = 0.
        self.x_const_arr = []
        self.x_xp = 0.
        self.x_xp_err = 0.
        self.x_xp_arr = []
        self.x_yp = 0.
        self.x_yp_err = 0.
        self.x_yp_arr = []
        self.x_xp2 = 0.
        self.x_xp2_err = 0.
        self.x_xp2_arr = []
        self.x_xpyp = 0.
        self.x_xpyp_err = 0.
        self.x_xpyp_arr = []
        self.x_yp2 = 0.
        self.x_yp2_err = 0.
        self.x_yp2_arr = []
        self.x_xp3 = 0.
        self.x_xp3_err = 0.
        self.x_xp3_arr = []
        self.x_xp2yp = 0.
        self.x_xp2yp_err = 0.
        self.x_xp2yp_arr = []
        self.x_xpyp2 = 0.
        self.x_xpyp2_err = 0.
        self.x_xpyp2_arr = []
        self.x_yp3 = 0.
        self.x_yp3_err = 0.
        self.x_yp3_arr = []

        self.y_const = 0.
        self.y_const_err = 0.
        self.y_const_arr = []
        self.y_xp = 0.
        self.y_xp_err = 0.
        self.y_xp_arr = []
        self.y_yp = 0.
        self.y_yp_err = 0.
        self.y_yp_arr = []
        self.y_xp2 = 0.
        self.y_xp2_err = 0.
        self.y_xp2_arr = []
        self.y_xpyp = 0.
        self.y_xpyp_err = 0.
        self.y_xpyp_arr = []
        self.y_yp2 = 0.
        self.y_yp2_err = 0.
        self.y_yp2_arr = []
        self.y_xp3 = 0.
        self.y_xp3_err = 0.
        self.y_xp3_arr = []
        self.y_xp2yp = 0.
        self.y_xp2yp_err = 0.
        self.y_xp2yp_arr = []
        self.y_xpyp2 = 0.
        self.y_xpyp2_err = 0.
        self.y_xpyp2_arr = []
        self.y_yp3 = 0.
        self.y_yp3_err = 0.
        self.y_yp3_arr = []

    ## figure out how many fields are in the file
    def getCount(self,fname,reset=True):
        ## set counter
        self.N_fields = 0
        ## encoding is not 'utf-8' for some reason
        with open(fname,'r',encoding='iso-8859-1') as f:
            ## read in everything
            buff = f.read()
            ## count how many 'Start's there are
            self.N_fields = len(re.findall("Start",buff))

    ## separate out values
    def parseEqn(self,eqn,isxy):
        const = float(eqn.split()[2])
        xp = float(eqn.split()[3].split('*')[0])
        yp = float(eqn.split()[4].split('*')[0])
        xp2 = float(eqn.split()[5].split('*')[0])
        xpyp = float(eqn.split()[6].split('*')[0])
        yp2 = float(eqn.split()[7].split('*')[0])
        xp3 = float(eqn.split()[8].split('*')[0])
        xp2yp = float(eqn.split()[9].split('*')[0])
        xpyp2 = float(eqn.split()[10].split('*')[0])
        yp3 = float(eqn.split()[11].split('*')[0])
        if(isxy == 'x'):
            self.x_const_arr.append(const)
            self.x_xp_arr.append(xp)
            self.x_yp_arr.append(yp)
            self.x_xp2_arr.append(xp2)
            self.x_xpyp_arr.append(xpyp)
            self.x_yp2_arr.append(yp2)
            self.x_xp3_arr.append(xp3)
            self.x_xp2yp_arr.append(xp2yp)
            self.x_xpyp2_arr.append(xpyp2)
            self.x_yp3_arr.append(yp3)
            return 1
        elif(isxy == 'y'):
            self.y_const_arr.append(const)
            self.y_xp_arr.append(xp)
            self.y_yp_arr.append(yp)
            self.y_xp2_arr.append(xp2)
            self.y_xpyp_arr.append(xpyp)
            self.y_yp2_arr.append(yp2)
            self.y_xp3_arr.append(xp3)
            self.y_xp2yp_arr.append(xp2yp)
            self.y_xpyp2_arr.append(xpyp2)
            self.y_yp3_arr.append(yp3)
            return 2
        else:
            return 0


    ## read in Astrometrica output
    def readAstrometrica(self,fname,reset=True):
        if(reset):
            self.clearVars()
        self.getCount(fname,reset)
        ## encoding is not 'utf-8' for some reason
        with open(fname,'r',encoding='iso-8859-1') as f:
            ## read the file
            all_lines = f.readlines()

            ## step through each line in the file
            for i in range(len(all_lines)):
                if("Start" in all_lines[i]):
                    offset = 0
                    ## find offset
                    for j in range(16,22):
                        if "Astrometry" in all_lines[i+j]:
                            offset = j+1
                    ## grab the number of stars used
                    self.N_stars.append(int(all_lines[i+offset].split()[0]))

                    ## information for the X solution
                    self.x_eqn = all_lines[i+offset+1]+all_lines[i+offset+2]+all_lines[i+offset+3]
                    ## print the equation if we want to
                    if(self._verbose_):
                        print(self.x_eqn)
                    ## store the values from the equation
                    self.parseEqn(self.x_eqn,'x')
                    ## print coefficients if we want to
                    if(self._verbose_):
                        print(self.x_const,self.x_xp,self.x_yp,self.x_xp2,self.x_xpyp,self.x_yp2,self.x_xp3,self.x_xp2yp,self.x_xpyp2,self.x_yp3)

                    ## information for the Y solution
                    self.y_eqn = all_lines[i+offset+4]+all_lines[i+offset+5]+all_lines[i+offset+6]
                    ## print the equation if we want to
                    if(self._verbose_):
                        print(self.y_eqn)
                    ## store the values from the equation
                    self.parseEqn(self.y_eqn,'y')
                    ## print coefficients if we want to
                    if(self._verbose_):
                        print(self.y_const,self.y_xp,self.y_yp,self.y_xp2,self.y_xpyp,self.y_yp2,self.y_xp3,self.y_xp2yp,self.y_xpyp2,self.y_yp3)

        ## calculate the average and stdev values
        self.x_const = statistics.mean(self.x_const_arr)
        self.x_xp = statistics.mean(self.x_xp_arr)
        self.x_yp = statistics.mean(self.x_yp_arr)
        self.x_xp2 = statistics.mean(self.x_xp2_arr)
        self.x_xpyp = statistics.mean(self.x_xpyp_arr)
        self.x_yp2 = statistics.mean(self.x_yp2_arr)
        self.x_xp3 = statistics.mean(self.x_xp3_arr)
        self.x_xp2yp = statistics.mean(self.x_xp2yp_arr)
        self.x_xpyp2 = statistics.mean(self.x_xpyp2_arr)
        self.x_yp3 = statistics.mean(self.x_yp3_arr)

        self.y_const = statistics.mean(self.y_const_arr)
        self.y_xp = statistics.mean(self.y_xp_arr)
        self.y_yp = statistics.mean(self.y_yp_arr)
        self.y_xp2 = statistics.mean(self.y_xp2_arr)
        self.y_xpyp = statistics.mean(self.y_xpyp_arr)
        self.y_yp2 = statistics.mean(self.y_yp2_arr)
        self.y_xp3 = statistics.mean(self.y_xp3_arr)
        self.y_xp2yp = statistics.mean(self.y_xp2yp_arr)
        self.y_xpyp2 = statistics.mean(self.y_xpyp2_arr)
        self.y_yp3 = statistics.mean(self.y_yp3_arr)

        if(len(self.N_stars) > 1):
            self.x_const_err = statistics.stdev(self.x_const_arr)
            self.x_xp_err = statistics.stdev(self.x_xp_arr)
            self.x_yp_err = statistics.stdev(self.x_yp_arr)
            self.x_xp2_err = statistics.stdev(self.x_xp2_arr)
            self.x_xpyp_err = statistics.stdev(self.x_xpyp_arr)
            self.x_yp2_err = statistics.stdev(self.x_yp2_arr)
            self.x_xp3_err = statistics.stdev(self.x_xp3_arr)
            self.x_xp2yp_err = statistics.stdev(self.x_xp2yp_arr)
            self.x_xpyp2_err = statistics.stdev(self.x_xpyp2_arr)
            self.x_yp3_err = statistics.stdev(self.x_yp3_arr)

            self.y_const_err = statistics.stdev(self.y_const_arr)
            self.y_xp_err = statistics.stdev(self.y_xp_arr)
            self.y_yp_err = statistics.stdev(self.y_yp_arr)
            self.y_xp2_err = statistics.stdev(self.y_xp2_arr)
            self.y_xpyp_err = statistics.stdev(self.y_xpyp_arr)
            self.y_yp2_err = statistics.stdev(self.y_yp2_arr)
            self.y_xp3_err = statistics.stdev(self.y_xp3_arr)
            self.y_xp2yp_err = statistics.stdev(self.y_xp2yp_arr)
            self.y_xpyp2_err = statistics.stdev(self.y_xpyp2_arr)
            self.y_yp3_err = statistics.stdev(self.y_yp3_arr)

    def printResults(self,fname,desc='',append=False):
        openType = 'a' if append else 'w'
        with open(fname,openType) as f:
            if(not append):
                f.write('desc,x_const,x_const_err,x_xp,x_xp_err,x_yp,x_yp_err,x_xp2,x_xp2_err,x_xpyp,x_xpyp_err,x_yp2,x_yp2_err,x_xp3,x_xp3_err,x_xp2yp,x_xp2yp_err,x_xpyp2,x_xpyp2_err,x_yp3,x_yp3_err,y_const,y_const_err,y_xp,y_xp_err,y_yp,y_yp_err,y_xp2,y_xp2_err,y_xpyp,y_xpyp_err,y_yp2,y_yp2_err,y_xp3,y_xp3_err,y_xp2yp,y_xp2yp_err,y_xpyp2,y_xpyp2_err,y_yp3,y_yp3_err\n')
            f.write('%s,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e,%13.6e\n'%(desc,self.x_const,self.x_const_err,self.x_xp,self.x_xp_err,self.x_yp,self.x_yp_err,self.x_xp2,self.x_xp2_err,self.x_xpyp,self.x_xpyp_err,self.x_yp2,self.x_yp2_err,self.x_xp3,self.x_xp3_err,self.x_xp2yp,self.x_xp2yp_err,self.x_xpyp2,self.x_xpyp2_err,self.x_yp3,self.x_yp3_err,self.y_const,self.y_const_err,self.y_xp,self.y_xp_err,self.y_yp,self.y_yp_err,self.y_xp2,self.y_xp2_err,self.y_xpyp,self.y_xpyp_err,self.y_yp2,self.y_yp2_err,self.y_xp3,self.y_xp3_err,self.y_xp2yp,self.y_xp2yp_err,self.y_xpyp2,self.y_xpyp2_err,self.y_yp3,self.y_yp3_err))

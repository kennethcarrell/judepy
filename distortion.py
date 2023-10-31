####
#### a class to read in distortion information from Astrometrica ouput
####
class Distortion:
    ## give variables some initial value
    N_stars = 0
    x_eqn = y_eqn = ''
    x_const = x_xp = x_yp = x_xp2 = x_xpyp = x_yp2 = x_xp3 = x_xp2yp = x_xpyp2 = x_yp3 = 0.
    y_const = y_xp = y_yp = y_xp2 = y_xpyp = y_yp2 = y_xp3 = y_xp2yp = y_xpyp2 = y_yp3 = 0.
    _verbose_ = False

    ## class constructor
    def __init__(self,fname='',vb=False):
        self._verbose_ = vb
        if(fname != ''):
            self.readAstrometrica(fname)

    ## read in Astrometrica output
    def readAstrometrica(self,fname):
        ## encoding is not 'utf-8' for some reason
        with open(fname,'r',encoding='iso-8859-1') as f:
            ## read the file
            all_lines = f.readlines()
            ## grab the number of stars used
            self.N_stars = int(all_lines[19].split()[0])

            ## information for the X solution
            self.x_eqn = all_lines[20]+all_lines[21]+all_lines[22]
            ## print the equation if we want to
            if(self._verbose_):
                print(self.x_eqn)
            self.x_const = float(self.x_eqn.split()[2])
            self.x_xp = float(self.x_eqn.split()[3].split('*')[0])
            self.x_yp = float(self.x_eqn.split()[4].split('*')[0])
            self.x_xp2 = float(self.x_eqn.split()[5].split('*')[0])
            self.x_xpyp = float(self.x_eqn.split()[6].split('*')[0])
            self.x_yp2 = float(self.x_eqn.split()[7].split('*')[0])
            self.x_xp3 = float(self.x_eqn.split()[8].split('*')[0])
            self.x_xp2yp = float(self.x_eqn.split()[9].split('*')[0])
            self.x_xpyp2 = float(self.x_eqn.split()[10].split('*')[0])
            self.x_yp3 = float(self.x_eqn.split()[11].split('*')[0])
            ## print coefficients if we want to
            if(self._verbose_):
                print(self.x_const,self.x_xp,self.x_yp,self.x_xp2,self.x_xpyp,self.x_yp2,self.x_xp3,self.x_xp2yp,self.x_xpyp2,self.x_yp3)

            ## information for the Y solution
            self.y_eqn = all_lines[23]+all_lines[24]+all_lines[25]
            ## print the equation if we want to
            if(self._verbose_):
                print(self.y_eqn)
            self.y_const = float(self.y_eqn.split()[2])
            self.y_xp = float(self.y_eqn.split()[3].split('*')[0])
            self.y_yp = float(self.y_eqn.split()[4].split('*')[0])
            self.y_xp2 = float(self.y_eqn.split()[5].split('*')[0])
            self.y_xpyp = float(self.y_eqn.split()[6].split('*')[0])
            self.y_yp2 = float(self.y_eqn.split()[7].split('*')[0])
            self.y_xp3 = float(self.y_eqn.split()[8].split('*')[0])
            self.y_xp2yp = float(self.y_eqn.split()[9].split('*')[0])
            self.y_xpyp2 = float(self.y_eqn.split()[10].split('*')[0])
            self.y_yp3 = float(self.y_eqn.split()[11].split('*')[0])
            ## print coefficients if we want to
            if(self._verbose_):
                print(self.y_const,self.y_xp,self.y_yp,self.y_xp2,self.y_xpyp,self.y_yp2,self.y_xp3,self.y_xp2yp,self.y_xpyp2,self.y_yp3)

####
#### example usage:
####
#### f1 = Distortion('Bruns HIP 100453 ED127W.log')
#### f2 = Distortion('Bruns HIP 102488 RC8.log')
#### print(f1.N_stars,f1.x_const,f1.y_const)
#### print(f2.N_stars,f2.x_const,f2.y_const)
#### f3 = Distortion()
#### f3.readAstrometrica('Bruns HIP 100453 ED127W.log')
#### print(f3.N_stars,f3.x_const,f3.y_const)

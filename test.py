from distortion import Distortion
import tkinter as tk
from tkinter import filedialog
from tkinter.messagebox import askyesno
from os.path import basename

## create a fake main window
root = tk.Tk()
root.withdraw()

## ask for the filenames
fname = filedialog.askopenfilenames()

combine = False
if(len(fname)>1):
    combine = askyesno(title='Combine?', message='Would you like to combine all results from all files?')

## run each file
for name in fname:
    f1 = Distortion(name)
    print(f1.N_fields,f1.x_const,f1.x_const_err,f1.y_const,f1.y_const_err)
    append = True if len(f1.N_stars)>1 else False
    f1.printResults('output.csv',desc=basename(name),append=append)

if(combine):
    for name in fname:
        fcomb = Distortion(name,reset=False)
        desc = desc+'+'+basename(name) if len(fcomb.N_stars)>1 else basename(name)
        print(fcomb.N_fields,len(fcomb.N_stars),fcomb.x_const,fcomb.x_const_err,fcomb.y_const,fcomb.y_const_err)
        fcomb.printResults('output.csv',desc=desc,append=True)

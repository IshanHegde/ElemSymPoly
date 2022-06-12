#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 05:44:13 2022

@author: ishanhegde
"""

import tkinter as tk
import pandas as pd
import openpyxl
from tkinter.filedialog import askopenfilename, asksaveasfilename
from tkinter import *
import numpy as np
import itertools
import tensorflow as tf

class App(object):
    
    def __init__(self, master, **kwargs):
        
        self.master=master
        self.create_window()
        self.filepath=''
        self.text=''
        self.dataFrame=pd.DataFrame()
        self.itemVar=pd.DataFrame()
        self.personVar=pd.DataFrame()
        self.outputpath=''
    
    def create_window(self):
        
        self.master.title("Rasch Model")
        self.master.rowconfigure(0, minsize=800, weight=1)
        self.master.columnconfigure(1, minsize=800, weight=1)
        
        fr_buttons = tk.Frame(self.master, relief=tk.RAISED, bd=3)
        
        self.txt_edit = tk.Text(self.master)
        self.btn_open = tk.Button(fr_buttons, text="Open", command=self.open_file)
        self.btn_save = tk.Button(fr_buttons, text="Save As...", command=self.save_file)
        self.btn_read = tk.Button(fr_buttons, text="Process...", command=self.read_file)
        
        self.btn_open.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        self.btn_save.grid(row=2, column=0, sticky="ew", padx=5)
        self.btn_read.grid(row=1,column=0,sticky="ew",padx=5,pady=5)

        fr_buttons.grid(row=0, column=0, sticky="ns")
        self.txt_edit.grid(row=0, column=1, sticky="nsew")

    def open_file(self):
        
        self.filepath = askopenfilename(filetypes=[("Excel Files",".xlsx"),('CSV Files',"*.csv"),("Text Files", "*.txt"),("All Files", "*.*")])
        if not self.filepath:
            return 
        
        self.txt_edit.delete(1.0,tk.END)
        
        if self.filepath[-4:]=='.csv':
            self.dataFrame=pd.read_csv(self.filepath)
        elif self.filepath[-4:]=='.txt':
            self.dataFrame=pd.read_table(self.filepath)
        elif self.filepath[-5:]=='.xlsx':
            self.dataFrame=pd.read_excel(self.filepath,engine='openpyxl')
        else:
            self.txt_edit.insert(tk.END, "Unrecognized file format")
            return
        
        self.txt_edit.insert(tk.END,f"File Path: {self.filepath}")
        self.txt_edit.insert(tk.END,"\n"*10)
        self.txt_edit.insert(tk.END,"Item Variables\nPerson Variables\n")
        
        for i in range(len(self.dataFrame.columns)):
            self.txt_edit.insert(tk.END,f"{i+1}. {self.dataFrame.columns[i]}\n")
    
    
        self.master.title(f"Rasch Model- {self.filepath}")
        
    def read_file(self):
        
        self.text = self.txt_edit.get(1.0,tk.END).split('\n')
        
        i=0
        item_var_names=[]
        person_var_names=[]
        
        while i < len(self.text)-1:
            
            while i<len(self.text)-1 and self.text[i] !="Item Variables":
                i+=1
                
            i+=1
            
            while i<len(self.text)-1 and self.text[i] !="Person Variables":
                if self.text[i]!="" :
                    item_var_names.append(self.text[i][3:])
                i+=1
                
            i+=1
            
            while i<len(self.text)-1:
                if self.text[i]!="" :
                    person_var_names.append(self.text[i][3:])
                i+=1
        
        
        for var in item_var_names:
            self.itemVar = pd.concat([self.itemVar,self.dataFrame[var]],axis=1)
        
        for var in person_var_names:
            self.personVar = pd.concat([self.personVar,self.dataFrame[var]],axis=1)
        
        print(self.itemVar)
        print(self.personVar)
    
    def save_file(self):
        
        self.outputpath=asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt"),('CSV Files',"*.csv"),("All Files", "*.*")],)
        if not self.outputpath:
            return
        
        with open(self.outputpath, "w") as output_file:
            text = self.txt_edit.get(1.0, tk.END)
            output_file.write(text)
            self.master.title(f"Rasch Model - {self.outputpath}")
        
    

def sigmoid(x):
    return 1/(1+np.exp(-x))



def main():
    
    root = tk.Tk()
    app = App(root)
    root.mainloop()

if __name__ == "__main__":
    
    main()


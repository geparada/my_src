#!/usr/bin/python
# -*- coding: utf-8 -*-

#
##
###
### Script: multifasta2fasta.py
### Author: Jonathan Grandaubert
###
### Purpose: 
### Create fasta files from a multifasta file.
###
### Copyright (C) 2008 - Jonathan Grandaubert
### 
### It is a free software; you can redistribute it and/or
### modify it under the terms of the GNU General Public License
### as published by the Free Software Foundation; either version 2
### of the License, or (at your option) any later version.
### 
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###
##
#


from Tkinter import *
from tkFileDialog import askopenfilename,askdirectory
import re,string,os


def multifasta2fasta(file):
	try:
		fi=open(file,'r')
		content=fi.readlines()
		fi.close()

		for i in content:
			if re.search('>',i):
				l=re.sub(">|\n|\r","",i)
				l=string.split(l)[0]
				l=re.sub(",|;|:|!","",l)
				filename="%s/%s.fasta" %(re.sub('\n','',text_rep2.get(1.0,END)),l)
				fo=open(filename,'w')
			fo.write(i)
		fo.close()
	except IOError,OSError:
		pass

def openfile():
	Formats=[('All Files','*.*'),('Fasta Files','*.fasta'),('Text Files','*.txt')]
	filename=askopenfilename(parent=root,filetypes=Formats,title="Open File")
	if filename:
		text_rep.delete(1.0,END)
		text_rep.insert(END,filename)
		text_run.delete(1.0,END)
		

def opendirectory():
	directory=askdirectory(parent=root)
	if directory:
		text_rep2.delete(1.0,END)
		text_rep2.insert(END,directory)
		text_run.delete(1.0,END)		

def run():
	text_run.delete(1.0,END)
	filename=text_rep.get(1.0,END)
	filename=re.sub('\n','',filename)

	multifasta2fasta(filename)
	
	text_run.insert(END,"Done!")


root=Tk()
root.title("multifasta2fasta")
root.resizable(False,False)

frame_all=Frame(root)
frame_all.grid()

Label(frame_all).grid(row=0)
Label(frame_all,text="Input: multifasta file",fg="white",bg="black").grid(row=1,sticky=W+E,padx=5)
frame_input=Frame(frame_all,bd=3,relief=GROOVE)
frame_input.grid(row=2,column=0,sticky=W+E,padx=5)

frame_input_data=Frame(frame_input)
frame_input_data.grid(row=1)

bouton_rep=Button(frame_input_data,text="Open multifasta file",height=1,width=15,command=openfile)
bouton_rep.grid(row=1,column=0,sticky='W',padx=5,pady=5)
text_rep=Text(frame_input_data,bg="white",height=1,width=50)
text_rep.grid(row=1,column=1,sticky='W',padx=5,pady=5)

Label(frame_all).grid(row=3)
label_step2=Label(frame_all,text="Output: directory for fasta files",fg="white",bg="black").grid(row=4,sticky=W+E,padx=5)
frame_output=Frame(frame_all,bd=3,relief=GROOVE)
frame_output.grid(row=5,column=0,sticky=W+E,padx=5)

frame_output_data=Frame(frame_output)
frame_output_data.grid(row=1)

bouton_rep2=Button(frame_output_data,text="Save fasta files",height=1,width=15,command=opendirectory)
bouton_rep2.grid(row=1,column=0,sticky='W',padx=5,pady=5)
text_rep2=Text(frame_output_data,bg="white",height=1,width=50)
text_rep2.grid(row=1,column=1,sticky='W',padx=5,pady=5)

Label(frame_all).grid(row=6)
bouton_run=Button(frame_all,text="Convert",bg="grey75",height=1,width=15,command=run)
bouton_run.grid(row=7,column=0,sticky=E,padx=5)
text_run=Text(frame_all,height=1,width=10,bg="grey85",relief=FLAT)
text_run.grid(row=7,sticky=W,padx=5)
Label(frame_all).grid(row=8)


root.mainloop()

import tkFileDialog
import Pmw

import pymol
from pymol.cgo import * 
import os
import sys
import re

#------------------------------------------------------------------------------#
# Create Pharao Menu                                                           #
#------------------------------------------------------------------------------#
def __init__(self):
	self.menuBar.addmenuitem('Plugin', 'command',
							'Pharao visualization plugin', label = 'Pharao...',
							command = lambda s=self: PharaoPlugin(s))
											

class PharaoPlugin:

	#-----------------------
	def __init__(self, app):

		print "+++ SILICOS::PHARAO Pymol Plugin (2010) +++"
		self.phar = []
		#TODO store id for pharmacophore

		self.parent = app.root
		self.dialog = Pmw.Dialog(self.parent, \
								buttons = ( \
								'Create pharmacophore', 
								'Read pharmacophore...', 
								'Write pharmacophore...', 
								'Create exclusion spheres',
								'Quit'),
								title = 'Pharao Pymol plugin',
								buttonboxpos = 'e',
								command = self.execute)

	#----------------------
	def execute(self, opt):
		if(opt=='Create pharmacophore'):
			self.runPharao()
		elif(opt=='Read pharmacophore...'):
			self.readPharao()
		elif(opt=='Write pharmacophore...'):
			self.writePharao()
		elif(opt=='Create exclusion spheres'):
			self.createExcl()
		else:
			print "+++ SILICOS::PHARAO Pymol Plugin (2010) +++"
			self.dialog.withdraw()

	#------------------
	def log(self, msg):
		#TODO create label in widget
		print " * "+msg

	#TODO method 'err' ?

#------------------------------------------------------------------------------#
# Create Pharmacophore                                                         #
#------------------------------------------------------------------------------#

	def runPharao(self):

		#check if selection is made
		if(not("sele" in pymol.cmd.get_names("selections"))):
			self.log("ERROR: no selection 'sele' defined.")
			sys.exit(1)

		#get object-id from selection
		list = pymol.cmd.index("sele")
		id = list[0][0]
		for x in list:
			if x[0]!=id:
				self.log("ERROR: multiple selection!")
				sys.exit(1)
		
		#add hydrogens and save in SDF format
		pymol.cmd.h_add("sele")
		pymol.cmd.save("/tmp/_tmp_Pharao.pdb", id) 
		pymol.cmd.remove("hydrogens in "+id)
  
		#run Pharao
		os.system("/usr/local/bin/pharao -d /tmp/_tmp_Pharao.pdb -p /tmp/_tmp_Pharao.phar")
		file = open("/tmp/_tmp_Pharao.phar")
  	
		#display results  
		if file != None:
			self.parseFile(file)
			file.close()
		else:
			self.log("ERROR: can't read Pharmacophore file.")
		
		# remove temporary files 
		os.system("rm /tmp/_tmp_Pharao.pdb")
		os.system("rm /tmp/_tmp_Pharao.phar")

		self.log("Done creating pharmacophore.")

		pymol.cmd.center()
		pymol.cmd.zoom()
	
#------------------------------------------------------------------------------#
# Read Pharmacophore    																											 #
#------------------------------------------------------------------------------#

	def readPharao(self):

		file = tkFileDialog.askopenfile(parent=self.parent,
                                    mode='r',
									filetypes=[('Pharmacophore', '*.phar')],
                                    title='Open pharmacophore file')
		if file != None:
			self.parseFile(file)
			file.close()
		else:
			self.log("ERROR: can't read pharmacophore file.")

		self.log("Done reading pharmacophore.")
		
		pymol.cmd.center()
		pymol.cmd.zoom()
	
#------------------------------------------------------------------------------#
# Create Exclusion Spheres                                                     #
#------------------------------------------------------------------------------#

	def createExcl(self):

		#check if selection is made
		if(not("sele" in pymol.cmd.get_names("selections"))):
			self.log("ERROR: no selection 'sele' defined.")
			sys.exit(1)

		pymol.cmd.select("_tmp_1_", "sele around 5 & !sele")

		exclPoint = [ COLOR, 0.25, 0.25, 0.25 ]

		c = 0
		model = pymol.cmd.get_model("_tmp_1_")
		for a in model.atom:
			# add point to current phar
			point = ("EXCL", a.coord[0], a.coord[1], a.coord[2], 1.4, "0", 0.0, 0.0, 0.0)
			self.phar.append(point)
			
			exclPoint.extend([ SPHERE, a.coord[0], a.coord[1], a.coord[2], 1.4])
			c = c + 1

		if(len(exclPoint) != 0 ):
			pymol.cmd.load_cgo(exclPoint,"excl",1)
			pymol.cmd.set('cgo_transparency',0.2,"excl")
		
		pymol.cmd.delete("_tmp_1_")
		self.log("Done Creating "+str(c)+" exclusion spheres.")
	
		pymol.cmd.center()
		pymol.cmd.zoom()

#------------------------------------------------------------------------------#
# Write Pharmacophore                                                          #
#------------------------------------------------------------------------------#

	def writePharao(self):

		file = tkFileDialog.asksaveasfile(parent=self.parent,
										filetypes=[('Pharmacophore', '*.phar')],
										title='Save pharmacophore')
		if file == None:
			self.log("ERROR: cannot save file.")
			sys.exit(1)
		
		file.write('Generated_by_Pharao_Pymol_Plugin\n')
		for p in self.phar:
			if(p[5] != "0"):
				file.write(p[0]+'\t'+str(p[1])+'\t'+str(p[2])+'\t'+str(p[3])+'\t'+str(p[4])+'\t'+
									 str(p[5])+'\t'+str(p[6])+'\t'+str(p[7])+'\t'+str(p[8])+'\n')
			else:
				file.write(p[0]+'\t'+str(p[1])+'\t'+str(p[2])+'\t'+str(p[3])+'\t'+str(p[4])+'\t'+		
				                     str(p[5])+'\t'+"0.000"+'\t'+"0.000"+'\t'+"0.000"+'\n')
		file.write('$$$$\n')
		file.close()
		
		self.log("Done writing pharmacophore.")	

#------------------------------------------------------------------------------#
# parse .phar file                                                             #
#------------------------------------------------------------------------------#

	def parseFile(self, file):

		# set colors
		aromPoint = [ COLOR, 0.50, 0.85, 0.05 ]
		lipoPoint = [ COLOR, 0.50, 0.10, 0.05 ]
		hyblPoint = [ COLOR, 0.50, 0.50, 0.05 ]
		hdonPoint = [ COLOR, 0.20, 0.65, 0.85 ]
		haccPoint = [ COLOR, 0.80, 0.65, 0.85 ]
		hybhPoint = [ COLOR, 0.50, 0.65, 0.85 ]
		poscPoint = [ COLOR, 0.00, 0.00, 1.00 ]
		negcPoint = [ COLOR, 1.00, 0.00, 0.00 ]
		exclPoint = [ COLOR, 0.25, 0.25, 0.25 ]
		normal = []

		id = file.readline()
		id = id.rstrip()

		#read pharmacophore points
		for line in file:
			if re.match("\$",line) == None:
				strlist = line.split()
				if len(strlist) >= 9:
					x = float(strlist[1])
					y = float(strlist[2])
					z = float(strlist[3])
					sigma = 1.0/float(strlist[4])
					hasNormal = strlist[5]
					if(re.match("AROM",strlist[0]) != None):
						aromPoint.extend([SPHERE, x, y, z, sigma])

					elif(re.match("LIPO",strlist[0]) != None ):
						lipoPoint.extend([SPHERE, x, y, z, sigma])

					elif(re.match("HYBL",strlist[0]) != None):
						hyblPoint.extend([SPHERE, x, y, z, sigma])

					elif(re.match("HDON",strlist[0]) != None):
						hdonPoint.extend([SPHERE, x, y, z, sigma])
						
					elif(re.match("HACC",strlist[0]) != None):
						haccPoint.extend([SPHERE, x, y, z, sigma])
            
					elif(re.match("HYBH",strlist[0]) != None):
						hybhPoint.extend([SPHERE, x, y, z, sigma])
          
					elif(re.match("POSC",strlist[0]) != None ):
						poscPoint.extend([SPHERE, x, y, z, sigma])
        	
					elif(re.match("NEGC",strlist[0]) != None ):
						negcPoint.extend([SPHERE, x, y, z, sigma])
          
					elif(re.match("EXCL",strlist[0]) != None ):
						exclPoint.extend([SPHERE, x, y, z, sigma])

					else:
						self.log("WARNING: unknown point.")

					if(hasNormal == "1"):
						normal.extend( [ CYLINDER,
										x, y, z,
										float(strlist[6]), float(strlist[7]), float(strlist[8]),
										0.1,
										1.0, 1.0, 0.0,
										1.0, 1.0, 0.0
										])
						point = (strlist[0], x, y, z, sigma, hasNormal, float(strlist[6]), float(strlist[7]), float(strlist[8]))
						self.phar.append(point)
					else:
						point= (strlist[0], x, y, z, sigma, hasNormal, 0.0, 0.0, 0.0)
						self.phar.append(point)
						 
			else:
				self.log("pharmacophore file parsed correctly.")

		#load graphical objects into pymol 
		if ( len (normal) != 0 ):
			pymol.cmd.load_cgo(normal,"normals_"+id,1)
			pymol.cmd.set('cgo_transparency',0.0,"normals_"+id)
		if ( len (aromPoint) >= 5 ):
			pymol.cmd.load_cgo(aromPoint, "arom_"+id,1)
			pymol.cmd.set('cgo_transparency',0.4,"arom_"+id)
		if ( len (lipoPoint) >=5 ):
			pymol.cmd.load_cgo(lipoPoint, "lipo_"+id,1)
			pymol.cmd.set('cgo_transparency',0.4,"lipo_"+id)
		if ( len (hyblPoint) >=5 ):
			pymol.cmd.load_cgo(hyblPoint, "hybl_"+id,1)
			pymol.cmd.set('cgo_transparency',0.4,"hybl_"+id)
		if ( len (hdonPoint) >= 5 ):
			pymol.cmd.load_cgo(hdonPoint, "hdon_"+id,1)
			pymol.cmd.set('cgo_transparency',0.4,"hdon_"+id)
		if ( len (haccPoint) >= 5 ):
			pymol.cmd.load_cgo(haccPoint, "hacc_"+id,1)
			pymol.cmd.set('cgo_transparency',0.4,"hacc_"+id)
		if ( len (hybhPoint) >= 5 ):
			pymol.cmd.load_cgo(hybhPoint, "hybh_"+id,1)
			pymol.cmd.set('cgo_transparency',0.4,"hybh_"+id)
		if ( len (poscPoint) >= 5 ):
			pymol.cmd.load_cgo(poscPoint, "posc_"+id,1)
			pymol.cmd.set('cgo_transparency',0.4,"posc_"+id)
		if ( len (negcPoint) >= 5 ):
			pymol.cmd.load_cgo(negcPoint, "negc_"+id,1)
			pymol.cmd.set('cgo_transparency',0.4,"negc_"+id)
		if ( len (exclPoint) >= 5 ):
			pymol.cmd.load_cgo(exclPoint, "excl_"+id,1)
			pymol.cmd.set('cgo_transparency',0.2,"excl_"+id)

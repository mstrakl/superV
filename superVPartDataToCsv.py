import glob, os
import csv

resdir="parares"

mydir=os.getcwd()
os.chdir(resdir)

csvData=[]
for fname in glob.glob("*"):

	t = fname.split('.')
	t = t[-1]
	t = float(t)

	readX=False
	readU=False
	readAge=False

	datFArr=[]
	datFArr.append(t)
	with open(fname) as f:

		line = f.readline()
		line = f.readline()
		datArr=line.split(',')

		cnt=0
		for dat in datArr:
			datFArr.append(float(dat))
			if (cnt == 5):
				break;
			cnt += 1

	csvData.append(datFArr)
		
os.chdir(mydir)

csvData = sorted(csvData,key=lambda x: x[0])


# name of csv file  
filename = "partData.csv"
fields = ['tstep', 'x', 'y', 'z', 'Ux', 'Uy', 'Uz']  
    
# writing to csv file  
with open(filename, 'w') as csvfile:  
    # creating a csv writer object  
    csvwriter = csv.writer(csvfile)  
        
    # writing the fields  
    csvwriter.writerow(fields)  
        
    # writing the data rows  
    csvwriter.writerows(csvData) 
		

	
	


	

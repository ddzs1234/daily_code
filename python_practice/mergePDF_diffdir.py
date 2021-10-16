
import os
from PyPDF2 import PdfFileMerger,PdfFileReader

rootDirectory=r"/media/ashley/zsash/linux/junior/sed/AGNFITTER_OLDVERSION_TXT/AGNfitter-master/OUTPUT/"

foldersToSearch=[x[0] for x in os.walk(rootDirectory)]

pdfs=[]

for f in foldersToSearch:
    pdfs+=[x for x in os.listdir(f) if x.startswith("SED_manyrealizations")]
    
pdfs=list(set(pdfs))


# create the dictionary with each PDF file name as the keys and an empty list as the value#
pdfDict = dict((pdf, []) for pdf in pdfs)

# add the full paths to the pdfs that match the key name to each value list
for f in foldersToSearch:
    for item in os.listdir(f):
        if item.startswith("SED_manyrealizations"):
#            print 'item,pdf',item
#            print os.path.join(f,item)
            pdfDict[item]=(os.path.join(f,item))
#print pdfDict
       
merge=PdfFileMerger()
for item in pdfDict.keys():
    
    merge.append(PdfFileReader(file(pdfDict[item],'rb')))
merge.write('./all.pdf')
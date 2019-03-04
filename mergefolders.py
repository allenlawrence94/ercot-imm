import os
import shutil
import sys
publicdata = sys.argv[1]
securedata = sys.argv[2]
for item in os.walk(securedata):
    for filename in item[2]:
        filedate = filename.split('.')[2]
        newfile = os.path.join(publicdata,filedate,filename)
        oldfile = os.path.join(securedata,filename)
        shutil.move(oldfile,newfile)

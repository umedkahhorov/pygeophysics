import os

path = os.getcwd()
filenames = os.listdir(path)

for filename in filenames:
    os.rename(filename, filename.replace('kr0917_A3obs_s','').replace('-1.sgy','').replace('-00','-').replace('-0','-'))

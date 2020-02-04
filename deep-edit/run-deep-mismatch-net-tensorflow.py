"""

Usage
python3 run-deep-mismatch-net-tensorflow.py local_folder high_res

"""

import sys
import numpy as np
from scipy.sparse import coo_matrix
import glob
from keras import models, layers
import keras
from scipy import signal
from keras import backend as K
import os
import tensorflow as tf
import random
#from matplotlib.colors import LinearSegmentedColormap
#from matplotlib import pyplot as plt
#import matplotlib.pyplot as plt 
#import seaborn as sns
from keras.models import load_model


#%tensorflow_version 1.x
import tensorflow as tf
print(tf.__version__)

model = load_model(sys.argv[3])

# todo make param
resolution = int(sys.argv[2])
local_folder = sys.argv[1]

newLength = 500
colorscale = newLength

outputIs1Less = False

def clipAndScale(data, scalar):
  threshold = np.max(np.max(data))/scalar
  data[data > threshold] = threshold
  return data/threshold

def clipByDiagAndScale(data, scalar):
  threshold = np.mean(np.diagonal(data))/scalar
  data[data > threshold] = threshold
  return data/threshold

def getThresholdedData(matrix, scalar, diagtype):
  if diagtype:
    return clipByDiagAndScale(matrix, scalar)
  else:
    return clipAndScale(matrix, scalar)

def getDataThresholded(foldername, filepath, scalar, diagtype = True):
  matrix = np.load(foldername+'/'+filepath+'.npy')
  threshData = getThresholdedData(matrix, scalar, diagtype)
  return threshData

results = set()

def stringIsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def readInDataFromFoldersAndPredictBreaks(zipFileLocation, matrixFileLocations):
  with open(zipFileLocation) as f:
    content = [x.strip() for x in f.readlines()]
    numFiles = len(content)
    numFilesPerIter = 500
    numTimesToIter = int(numFiles/numFilesPerIter) + 1
    for runIter in range(numTimesToIter):
      startIndex = runIter * numFilesPerIter
      endIndex = min(startIndex+numFilesPerIter, numFiles)
      numFilesThisRun = endIndex - startIndex
      finData   = np.zeros([numFilesThisRun, newLength, newLength, 1])
      finIndices = np.zeros([numFilesThisRun, 2, 1, 1])
      for k in range(numFilesThisRun):
        fileK = content[k+startIndex]
        splitName = fileK.split('_')
        if stringIsInt(splitName[2]):
          finIndices[k,0,0,0] = int(splitName[2])
          finIndices[k,1,0,0] = int(splitName[4])
        else:
          finIndices[k,0,0,0] = int(splitName[3])
          finIndices[k,1,0,0] = int(splitName[5])
        finData[k, :, :, 0] = getDataThresholded(matrixFileLocations, fileK, colorscale)
        
      prediction = np.around(model.predict(finData))
      # ignore the last zero
      prediction = prediction[:,0:-1,:,:]
      for k in range(numFilesThisRun):
        x = np.squeeze(prediction[k,:,0,0])
        contiguous = int(finIndices[k,1,0,0] - finIndices[k,0,0,0]) == int(newLength/2)
        if(contiguous):
          #print('was contig')
          if np.count_nonzero(x) > 0:
            for y in np.nditer(np.nonzero(x)):
              results.add(int(finIndices[k,0,0,0])+y)
              #print(int(finIndices[k,0,0,0]),y)
        else:
          midpt = int((newLength-1)/2)
          #print('was not contig',finIndices[k,:,0,0])
          x1 = np.squeeze(prediction[k,0:midpt,0,0])
          x2 = np.squeeze(prediction[k,midpt+1:,0,0])
          if np.count_nonzero(x1) > 0:
            for y in np.nditer(np.nonzero(x1)):
              results.add(int(finIndices[k,0,0,0])+y)
              #print(int(finIndices[k,0,0,0]),y)
          if np.count_nonzero(x2) > 0:
            for y in np.nditer(np.nonzero(x2)):
              results.add(int(finIndices[k,1,0,0])+y)
              #print(int(finIndices[k,1,0,0]),y)

namesFileLocation1 = local_folder+'/chr_assembly_'+str(resolution)+'.txt'
actualMatrixFileLocations1 = local_folder+'/neural_net_run_'+str(resolution)
namesFileLocation2 = local_folder+'/diag_chr_assembly_'+str(resolution)+'.txt'
actualMatrixFileLocations2 = local_folder+'/neural_net_diag_run_'+str(resolution)

readInDataFromFoldersAndPredictBreaks(namesFileLocation1, actualMatrixFileLocations1)
readInDataFromFoldersAndPredictBreaks(namesFileLocation2, actualMatrixFileLocations2)

outputBreakPointsFile = open("mismatch_narrow.bed", "w")
results = list(results)
results.sort()
bufferForIndivCall = int(int(resolution)/4)
for index in results:
  # write line to output file
  genomicPosition = resolution*(index+1)
  outputBreakPointsFile.write("assembly")
  outputBreakPointsFile.write("\t")
  outputBreakPointsFile.write(str(genomicPosition - bufferForIndivCall))
  outputBreakPointsFile.write("\t")
  outputBreakPointsFile.write(str(genomicPosition + bufferForIndivCall))
  outputBreakPointsFile.write("\n")
outputBreakPointsFile.close()

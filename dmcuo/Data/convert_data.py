"""
Note: This file is just for me to convert the binary .mat blobs to .json format
"""
import numpy as np
import codecs, json
import numpy as np
import scipy.io
import matplotlib.pyplot as plt

save_path = "./dissolution_profiles_cu2p.json"  ## your path variable

# Convert data:
mat = scipy.io.loadmat('./dissolution_profiles.mat')
dissolution_profiles = mat['dissolution_profiles']
dissolution_profiles = np.array(dissolution_profiles)
print(dissolution_profiles.shape)
print(dissolution_profiles[0][0].shape)
data_time = np.array(dissolution_profiles[0][0][0])
print(data_time)
print(data_time.shape)

# Plots to find which set is mean/min/max:
if 0:
    print(dissolution_profiles)
    print(dissolution_profiles.shape)
    print(dissolution_profiles[0].shape)
    print(dissolution_profiles[0][0].shape)
    print(dissolution_profiles[0][0][0])
    print(dissolution_profiles[0][0][1])
    print(dissolution_profiles[0][0][2])
    print(dissolution_profiles[0][0][3])
    ## Time in h

    for ind in range(1,4):
        x = dissolution_profiles[0][0][0]
        y = dissolution_profiles[0][0][ind]
        print(x.shape)
        print(y.shape)
        for i in range(4):
            plt.plot(x[i], y[i], '*-')

    x = dissolution_profiles[0][0][0]
    y = dissolution_profiles[0][0][3]  # Max
    y = dissolution_profiles[0][0][1]  # Mean
    y = dissolution_profiles[0][0][2]  # Min

data = {'time': dissolution_profiles[0][0][0].tolist(),
        'mean': dissolution_profiles[0][0][1].tolist(),
        'min': dissolution_profiles[0][0][2].tolist(),
        'max': dissolution_profiles[0][0][3].tolist(),
        }
#c_cu2p

save_path
json.dump(data, codecs.open(save_path, 'w', encoding='utf-8'),
          separators=(',', ':'),
          sort_keys=True,
          indent=4)  ### this saves the array in .json format
if 0:
    a = np.arange(10).reshape(2, 5)  # a 2 by 5 array
    b = a.tolist()  # nested lists with same data, indices
    file_path = "/path.json"  ## your path variable
    json.dump(b, codecs.open(file_path, 'w', encoding='utf-8'),
              separators=(',', ':'),
              sort_keys=True,
              indent=4)  ### this saves the array in .json format


    obj_text = codecs.open(file_path, 'r', encoding='utf-8').read()
    b_new = json.loads(obj_text)
    a_new = np.array(b_new)
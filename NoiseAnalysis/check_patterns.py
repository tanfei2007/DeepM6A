import numpy as np
import h5py

f = h5py.File('train.hdf5')
X = np.array(f['x_train'])
y = np.array(f['y_train'])
idx_pos = np.where(y==1)
X_pos = X[idx_pos]

count_right1 = 0
count_right2 = 0
count_right = 0
dicts = {0: "A", 1: "C", 2: "G", 3: "T"}

def is_mutually_exclusive(lst):
    seen = []
    for i in lst:
        if i in seen:
            return False
        seen.append(i)
    return True

for i in X_pos:
    set1 = i[(31,33,35,37),:]
    set2 = i[(23,25,27,29),:]

    oset1 = []
    oset2 = []

    for a in set1:
        ans = dicts[np.argmax(a)]
        oset1.append(ans)
    for b in set2:
        ans = dicts[np.argmax(b)]
        oset2.append(ans)
    print oset2
    oset1ans = is_mutually_exclusive(oset1)
    #oset2ans = is_mutually_exclusive(oset2)
    oset2ans = (oset2 == ['T', 'G', 'C', 'A'])
    if oset1ans:
        print("First")
	count_right1 += 1
    if oset2ans:
        print("Second")
	count_right2 += 1
    if oset1ans and oset2ans:
        count_right += 1
print("Accuracy1:" + str(count_right1))
print("Accuracy2:" + str(count_right2))
print("Accuracy:" + str(count_right))

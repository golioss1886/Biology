import os,glob
import numpy as np
import functools
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
# from tensorflow.keras.models import Sequential
# from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense

def dir_cmp(x, y):
    if int(x.split('/')[-1].split('.')[0]) > int(y.split('/')[-1].split('.')[0]):
        return 1
    elif int(x.split('/')[-1].split('.')[0]) < int(y.split('/')[-1].split('.')[0]):
        return -1
    return 0

def SVM(X_train, X_test, y_train, y_test):
    # clf = svm.SVC(kernel='linear') # You can change the kernel to 'poly', 'rbf', etc.
    clf = svm.SVC(kernel='poly', degree=4) # You can change the kernel to 'poly', 'rbf', etc.
    # clf = svm.SVC(kernel='rbf') # You can change the kernel to 'poly', 'rbf', etc.
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    print(classification_report(y_test, y_pred))

def RF(X_train, X_test, y_train, y_test):
    RFscaler = RandomForestClassifier()
    RFscaler.fit(X_train,y_train)
    print(classification_report(RFscaler.predict(X_test),y_test))

def CNN():
    X_buff=X_train.reshape(1, len(X_train), 1)
    # 建立模型
    CNNmodel = Sequential()
    CNNmodel.add(Conv1D(32, 3, activation='relu', input_shape=(len(x), 1)))
    CNNmodel.add(MaxPooling1D(2))
    CNNmodel.add(Flatten())
    CNNmodel.add(Dense(1, activation='linear'))

    # 編譯模型
    CNNmodel.compile(loss='mean_squared_error', optimizer='adam')
    CNNmodel.fit(X_train,y_train, epochs=10, batch_size=1)
    print(metrics.classification_report(CNNmodel.predict(X_test),y_test))


def main():
    # Training
    np_files = glob.glob(os.path.join("train_np", '*'))
    np_files.sort(key=functools.cmp_to_key(dir_cmp))
    # for i in np_files:
        # print(i)
    stack_list = []
    for np_file in np_files:
        with open(np_file, 'rb') as f:
            data = np.load(f)
            stack_list.append(data)
    X = np.array(stack_list)
    X = X.reshape((80, -1))
    # print("\n")
    # print(X[-1])
    # for i in X:
        # print(i)
    y = np.zeros(80)  # create an array of zeros of size 80
    y[50:] = 1  # set the last 30 elements to 1
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # SVM(X_train, X_test, y_train, y_test)
    RF(X_train, X_test, y_train, y_test)
    # CNN()

if __name__ == '__main__':
    main()

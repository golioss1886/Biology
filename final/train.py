import os,glob
import numpy as np
import functools
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Flatten, Dense, Dropout, ZeroPadding2D
from tensorflow.keras.utils import to_categorical
# from keras.models import Sequential
# from keras.layers import Conv2D, MaxPooling2D, Flatten, Dense
os.environ['CUDA_VISIBLE_DEVICES'] = '4'


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

def CNN(X_train, X_test, y_train, y_test):
    # X_buff=X_train.reshape(1, len(X_train), 1)
    X_train = np.expand_dims(X_train, -1)
    y_train = np.expand_dims(y_train, -1)
    y_train = to_categorical(y_train, num_classes=2)


    X_test = np.expand_dims(X_test, -1)
    y_test = np.expand_dims(y_test, -1)
    # print(X_train.shape)
    # print(X_train)
    print(y_train.shape)
    '''
    # 建立模型
    model = Sequential()
    model.add(Conv1D(32, 3, activation='relu', input_shape=(len(x), 1)))
    model.add(MaxPooling1D(2))
    model.add(Flatten())
    model.add(Dense(1, activation='linear'))

    # 編譯模型
    model.compile(loss='mean_squared_error', optimizer='adam')
    model.fit(X_train,y_train, epochs=10, batch_size=1)
    print(metrics.classification_report(CNNmodel.predict(X_test),y_test))
    '''

    model = Sequential()
    model.add(ZeroPadding2D(padding=(1, 1), input_shape=(20,20,1)))
    # model.add(Conv2D(32, (3, 3), activation='relu', input_shape=(20, 20, 1)))
    model.add(Conv2D(32, (3, 3), activation='relu', input_shape=(22,22,1)))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    
    # model.add(ZeroPadding2D(padding=(1, 1), input_shape=(20,20,1)))
    model.add(ZeroPadding2D(padding=(1, 1), input_shape=(12,12,1)))
    model.add(Conv2D(64, (3, 3), activation='relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))

    model.add(ZeroPadding2D(padding=(1, 1), input_shape=(7,7,1)))
    model.add(Conv2D(128, (3, 3), activation='relu'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    # Add a convolutional layer
    # model.add(Conv2D(32, (20, 20), activation='relu', input_shape=(20, 20, 1)))
    # Add a pooling layer
    # model.add(MaxPooling2D())
    # Add another convolutional layer
    # model.add(Conv2D(64, (10, 10), activation='relu'))
    # model.add(Conv2D(128, (3, 3), activation='relu'))
    # model.add(MaxPooling2D(pool_size=(2, 2)))
    # Add another pooling layer
    # model.add(MaxPooling2D(pool_size=(2, 2)))
    # Flatten the tensor output from the previous layer
    model.add(Flatten())
    # Add a fully connected layer
    # model.add(Dense(512, activation='relu'))
    model.add(Dropout(0.2))  # Dropout layer with dropout rate of 0.5
    model.add(Dense(512, activation='relu'))
    model.add(Dense(128, activation='relu'))
    # Add the output layer with 10 neurons (for the 10 digits)
    # model.add(Dense(2, activation='softmax')) #?
    model.add(Dense(2, activation='softmax'))
    # Compile the model
    model.compile(loss='categorical_crossentropy',
    # model.compile(loss='binary_crossentropy',
    # model.compile(loss='binary_crossentropy',
                  optimizer='adam',
                  metrics=['accuracy'])
    print(X_train.shape)
    print(y_train.shape)
    model.fit(X_train,y_train, epochs=40, batch_size=2)
    # score = model.evaluate(X_test, y_test, verbose=0)
    # print('\n', 'Test accuracy:', score[1])
    predictions = model.predict(X_test)
    print(predictions.shape, "::::::::")
    answer = np.argmax(predictions, axis=-1)
    print(answer)
    accuracy = np.mean(answer == y_test)
    print(accuracy)

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
    Z = np.array(stack_list)
    print(Z.shape)
    X = Z.reshape((80, -1))
    # print("\n")
    # print(X[-1])
    # for i in X:
        # print(i)
    y = np.zeros(80)  # create an array of zeros of size 80
    y[50:] = 1  # set the last 30 elements to 1
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # SVM(X_train, X_test, y_train, y_test)
    # RF(X_train, X_test, y_train, y_test)
    X_train, X_test, y_train, y_test = train_test_split(Z, y, test_size=0.2, random_state=42)
    print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
    CNN(X_train, X_test, y_train, y_test)

if __name__ == '__main__':
    main()

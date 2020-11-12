import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.models import Sequential
from tensorflow.keras import layers
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l2
from sklearn.utils import shuffle
from sklearn.model_selection import KFold

def seq2vec(string, outlen, vocab):
    vector = [vocab[amino_acid] for amino_acid in string]
    vector = np.pad(vector, (0,outlen-len(vector)), constant_values=20)
    return np.array(vector)

def label(row):
    if row['binding score'] > 0:
        if row['digest score'] > 0: return 0
        else: return 2      
    elif row['binding score'] <= 0:
        if row['digest score'] <= 0: return 1
        else: return 3

#load data
data = pd.read_csv("family_101F.csv")
vfile = pd.read_csv("vocab.csv", index_col=None, skiprows=0)
vocab = {vfile["One-Letter Code"][i]:i for i in range(len(vfile["One-Letter Code"]))}

#change vector to one-hot encoding
Tx = data["len"].max()
vector = np.array([to_categorical(seq2vec(data["seq"][i], Tx, vocab)) for i in data.index])

#label the class according to scores
data.loc[(data["binding score"] + data["digest score"] > 0), "label"] = 0
data.loc[(data["binding score"] + data["digest score"] <= 0), "label"] = 1

#set up constant
nbatch = 320
nepoch = 200
kfold = 5
adam = Adam(learning_rate=0.0001)

kf = KFold(n_splits=kfold, shuffle=True, random_state=1)
acc, val_acc, loss, val_loss = pd.DataFrame([]), pd.DataFrame([]), pd.DataFrame([]), pd.DataFrame([])

for (train_index, test_index), k in zip(kf.split(vector,data["label"]), range(kfold)):
    
    #create training/testing set with k-fold    
    trainX = vector[train_index]
    trainY = data["label"][train_index].values
    testX = vector[test_index]
    testY = data["label"][test_index].values

    #train the model
    model = Sequential()
    model.add(layers.Bidirectional(layers.SimpleRNN(units=128, kernel_regularizer=l2(0.001), recurrent_regularizer=l2(0.001), input_shape=(Tx, len(vocab)), return_sequences=True)))
    model.add(layers.Bidirectional(layers.SimpleRNN(units=128, kernel_regularizer=l2(0.001), recurrent_regularizer=l2(0.001), input_shape=(Tx, len(vocab)), return_sequences=True)))
    model.add(layers.Bidirectional(layers.SimpleRNN(units=128, kernel_regularizer=l2(0.001), recurrent_regularizer=l2(0.001), input_shape=(Tx, len(vocab)))))
    #model.add(layers.Dropout(0.1, input_shape=(128,)))
    model.add(layers.Dense(units=128, activation='relu'))
    #model.add(layers.Dropout(0.1, input_shape=(128,)))
    model.add(layers.Dense(units=128, activation='relu'))
    #model.add(layers.Dropout(0.1, input_shape=(128,)))
    model.add(layers.Dense(units=2, activation='softmax'))
    model.compile(loss='sparse_categorical_crossentropy', metrics=['accuracy'], optimizer=adam)
    history = model.fit(trainX, trainY, validation_data=(testX,testY), epochs=nepoch, batch_size=nbatch, verbose=1)
    
    #keep model weights and result
    acc["k"+str(k+1)] = history.history["accuracy"]
    loss["k"+str(k+1)] = history.history["loss"] 
    val_acc["k"+str(k+1)] = history.history["val_accuracy"] 
    val_loss["k"+str(k+1)] = history.history["val_loss"] 
    model.save('saved_model/model-k'+str(k)) 

    #clean session
    del model
    tf.keras.backend.clear_session()
    
acc.to_csv("acc.csv", index_label=False)
loss.to_csv("loss.csv", index_label=False)
val_acc.to_csv("val_acc.csv", index_label=False)
val_loss.to_csv("val_loss.csv", index_label=False)

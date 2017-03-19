
# coding: utf-8

# In[1]:
import numpy as np
import tensorflow as tf
from tensorflow.contrib.tensorboard.plugins import projector
import time
import os


# In[2]:

class DataManager:

    def __init__(self, pathToInput, pathToLabels):
        self.input = np.genfromtxt(pathToInput, dtype = float, delimiter = ',', skip_header = 1)
        self.labels =  np.genfromtxt(pathToLabels, dtype=int, delimiter=',', skip_header=1)
        assert(len(self.input) == len(self.labels))
        
    #1 of k encoder for labels
    def onehot(self, target):
        return np.eye(len(np.unique(target)))[target]

    #creates training- & test-dataset
    #trainingShare: percental size of training-dataset
    def createDatasets(self, trainingShare):
        p = np.random.permutation(len(self.input))
        shuffled_input = self.input[p]
        shuffled_labels = self.labels[p]
        shuffled_labels_oh = self.onehot(shuffled_labels)
        
        trainingSetSize = int(trainingShare*len(self.input))

        tr_X, te_X = shuffled_input[:trainingSetSize, :], shuffled_input[trainingSetSize:, :]
        tr_Y, te_Y = shuffled_labels_oh[:trainingSetSize, :], shuffled_labels_oh[trainingSetSize:,:]

        return tr_X, te_X, tr_Y, te_Y, shuffled_labels[trainingSetSize:]
    
    #creates training- & test-dataset keeping the ratio of 0 to 1 in training and test-set equal
    def createBinaryDatasets(self, trainingShare):
        index = np.argsort(self.labels)
        sorted_in = self.input[index]
        
        zeroes = np.bincount(self.labels)[0]
        negatives = sorted_in[:zeroes,:]
        positives = sorted_in[zeroes:,:]
        
        p = np.random.permutation(negatives.shape[0])
        n_tr_X, n_te_X = np.split(negatives[p], [int(negatives.shape[0]*0.6)+1])

        p = np.random.permutation(positives.shape[0])
        p_tr_X, p_te_X = np.split(positives[p], [int(positives.shape[0]*0.6)])
                
        tr_X = np.vstack((n_tr_X, p_tr_X))
        te_X = np.vstack((n_te_X, p_te_X))
                
        tr_Y = np.zeros(tr_X.shape[0], dtype=int)
        tr_Y[n_tr_X.shape[0]:] = 1
        
        te_Y  = np.zeros(te_X.shape[0], dtype=int)
        te_Y[n_te_X.shape[0]:] = 1
        
        p = np.random.permutation(tr_Y.shape[0])
        p2 = np.random.permutation(te_Y.shape[0])
                
        return tr_X[p], te_X[p2], self.onehot(tr_Y[p]), self.onehot(te_Y[p2]), te_Y[p2]
        


# In[3]:

class Model():
    
    def __init__(self, input, h_sizes, output, learning_rate=0.01, p_keep_input=1.0, p_keep_hidden=1.0,l2_loss=False):
        self.h_sizes = h_sizes
        self.X = input
        self.Y = output
        self.lr = learning_rate
        self.pki = p_keep_input
        self.pkh= p_keep_hidden
        self.globalStep = globalStep
        self.l2_loss = l2_loss

    def bake(self):
        #create model
        py_x, py_x_pred,reguls = self.createModel()
        with tf.name_scope("cost"):
            #define cost function
            cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=py_x, labels=self.Y))
            #initialise regularisers
            if self.l2_loss:
                for r in reguls:
                    cost = cost+0.01*r
            training_cost = tf.summary.scalar("TrainingCost", cost)
            self.test_cost = tf.summary.scalar("TestCost", cost)
        with tf.name_scope("GradientDescent"):
            #define backprop alogrithm
            train_op = tf.train.GradientDescentOptimizer(self.lr).minimize(cost, global_step=globalStep)

        with tf.name_scope("predict"):
            predict_op = tf.argmax(py_x_pred, 1)

        with tf.name_scope('accuracy'):
            acc = tf.equal(predict_op, tf.argmax(self.Y, 1))
            acc = tf.reduce_mean(tf.cast(acc, tf.float32))
            training_acc = tf.summary.scalar("TrainingAccuracy", acc)
            self.test_acc = tf.summary.scalar("TestAccuracy", acc)
        
        #return training operator, prediction operator and probabilities for predictions
        return train_op, predict_op, tf.nn.softmax(py_x_pred)

    def init_weights(self, shape):
        return tf.Variable(tf.random_normal(shape, stddev=0.01))

    def createModel(self):
        ###Build weights
        with tf.name_scope("WeightsInput"):
            w_i = self.init_weights([int(self.X.shape[1]), self.h_sizes[0]])
            self.variable_summaries(w_i)
        w_matrices = [w_i]
        regularizers = [tf.nn.l2_loss(w_i)]
        with tf.name_scope("BiasInput"):
            b_i = tf.Variable(tf.zeros([self.h_sizes[0]]))
            self.variable_summaries(b_i)
        b_vectors = [b_i]

        for i in range(1,len(self.h_sizes)):
            with tf.name_scope("WeightsHidden%s"%i):
                w = self.init_weights([int(w_matrices[i - 1].get_shape()[1]), self.h_sizes[i]])
                self.variable_summaries(w)
            w_matrices.append(w)
            regularizers.append(tf.nn.l2_loss(w))
            with tf.name_scope("BiasHidden%s"%i):
                    b = tf.Variable(tf.zeros([self.h_sizes[i]]))
                    self.variable_summaries(b)
            b_vectors.append(b)

        with tf.name_scope("WeightsOutput"):
            w_o = self.init_weights([int(w_matrices[-1].get_shape()[1]), int(self.Y.shape[1])])
            self.variable_summaries(w_o)
        w_matrices.append(w_o)
        regularizers.append(tf.nn.l2_loss(w_o))
        with tf.name_scope("BiasOutput"):
            b_o = tf.Variable(tf.zeros([int(self.Y.shape[1])]))
            self.variable_summaries(b_o)
        b_vectors.append(b_o)


        ###Build layers
        with tf.name_scope("InputLayer"):
            x_drop = tf.nn.dropout(self.X, self.pki)
            in_layer = tf.nn.relu(tf.matmul(x_drop, w_matrices[0])+b_vectors[0])
            self.variable_summaries(in_layer)
        h_layers = [in_layer]

        for i in range(1, len(self.h_sizes)):
            with tf.name_scope("HiddenLayer" + str(i)):
                h_drop = tf.nn.dropout(h_layers[i - 1], self.pkh)
                hl = tf.nn.relu(tf.matmul(h_drop, w_matrices[i])+b_vectors[i])
                self.variable_summaries(hl)
            h_layers.append(hl)

        with tf.name_scope("Logits"):
            out_drop = tf.nn.dropout(h_layers[-1], self.pkh)
            logits = tf.matmul(out_drop, w_matrices[-1])+b_vectors[-1]
            self.variable_summaries(logits)

        ###Build prediction model without dropout
        with tf.name_scope("PredictionModel"):
            in_layer_pred = tf.nn.relu(tf.matmul(X, w_matrices[0])+b_vectors[0])
            h_layers_pred = [in_layer_pred]

            for i in range(1, len(self.h_sizes)):
                hl = tf.nn.relu(tf.matmul(h_layers_pred[i - 1], w_matrices[i])+b_vectors[i])
                h_layers_pred.append(hl)

            logits_pred = tf.matmul(h_layers_pred[-1], w_matrices[-1])+b_vectors[-1]
        
        return logits, logits_pred, regularizers
    
    #summarise in tensorboard
    def variable_summaries(self, var):
        with tf.name_scope('summaries'):
            mean = tf.reduce_mean(var)
            tf.summary.scalar('mean', mean)
            with tf.name_scope('stddev'):
                stddev = tf.sqrt(tf.reduce_mean(tf.square(var - mean)))
            tf.summary.scalar('stddev', stddev)
            tf.summary.scalar('max', tf.reduce_max(var))
            tf.summary.scalar('min', tf.reduce_min(var))
            #tf.summary.histogram('histogram', var)


# In[4]:

dm = DataManager("./res/jc_SRP057500_python.csv", "./res/labels_SRP057500_python.csv")
total_tp = total_tn = total_fp = total_fn = 0
#runs
for z in range(100):
    tr_X, te_X, tr_Y, te_Y, raw_labels = dm.createBinaryDatasets(0.6)
    print "####################Run %d####################" % (z+1)

    size_x = tr_X.shape[1]
    size_y = tr_Y.shape[1]

    with tf.name_scope("Inputs"):
        X = tf.placeholder("float", [None,size_x])
        Y = tf.placeholder("float", [None, size_y])
        
    #learning rate decay, currently not used
    with tf.name_scope("LearningRate"):
        globalStep = tf.Variable(0)
        learning_rate = tf.train.exponential_decay(0.1, globalStep, 200, 0.95, staircase=True)
        tf.summary.scalar("learning_rate", learning_rate)
    with tf.name_scope("Embedding"):
        embedding = tf.Variable(te_X, trainable=False, name='embedding')

    #create model
    m = Model(X,[100,100],Y, learning_rate=0.0001, p_keep_input=.65, p_keep_hidden=.65, l2_loss=True)
    train_op, predict_op, yhat = m.bake()

    max_test_acc = 0
    best_pred = []
    with tf.Session() as sess:
        ###tensorboard
        merged = tf.summary.merge_all()

        log_dir = './log/run%s'%time.time()
        writer = tf.summary.FileWriter(log_dir, sess.graph)
        saver = tf.train.Saver([embedding])

        with open(os.path.join(log_dir, 'metadata.tsv'), 'w') as metadata:
            for row in raw_labels:
                metadata.write('%d\n' % row)
        metadata.close()

        config = projector.ProjectorConfig()
        embed= config.embeddings.add()
        embed.tensor_name = embedding.name
        embed.metadata_path = os.path.join(log_dir, 'metadata.tsv')
        projector.visualize_embeddings(writer, config)
        ###/tensorboard
        
        tf.global_variables_initializer().run()
        #epochs
        for epoch in range(300):
            for i in range(len(tr_X)):
                sess.run(train_op, feed_dict={X: tr_X[i: i + 1], Y: tr_Y[i: i + 1]})
            
            train_accuracy = np.mean(np.argmax(tr_Y, axis=1) ==
                                     sess.run(predict_op, feed_dict={X: tr_X, Y: tr_Y}))
            test_accuracy = np.mean(np.argmax(te_Y, axis=1) ==
                                    sess.run(predict_op, feed_dict={X: te_X, Y: te_Y}))
            print("Epoch = %d, train accuracy = %.2f%%, test accuracy = %.2f%%"
                  % (epoch + 1, 100. * train_accuracy, 100. * test_accuracy))
            #save best prediction of the run
            if(test_accuracy>max_test_acc):
                max_test_acc = test_accuracy
                best_pred = sess.run(yhat, feed_dict={X: te_X, Y: te_Y})
                
            ###tensorboard
            if (epoch % 10 == 0):
                saver.save(sess, os.path.join(log_dir, 'model.ckpt'), epoch)
            training_summary = sess.run(merged, feed_dict={X: tr_X, Y: tr_Y})
            writer.add_summary(training_summary, epoch)
            summ_test_acc, summ_test_cost = sess.run([m.test_acc, m.test_cost], feed_dict={X:te_X, Y: te_Y})
            writer.add_summary(summ_test_acc, epoch)
            writer.add_summary(summ_test_cost, epoch)
            ###/tensorboard
        writer.close()

    true_labels = np.argmax(te_Y, axis=1)
    model_pred = np.argmax(best_pred, axis=1)
    scores = np.around(best_pred[np.arange(model_pred.shape[0]), model_pred], decimals=2)
    #calculate true positives etc. for this run
    tp = np.sum(np.logical_and(model_pred == 1, true_labels == 1), dtype=float)
    tn = np.sum(np.logical_and(model_pred == 0, true_labels == 0))
    fp = np.sum(np.logical_and(model_pred == 1, true_labels == 0), dtype=float)
    fn = np.sum(np.logical_and(model_pred == 0, true_labels == 1))
    print 'TP: %i, FP: %i, TN: %i, FN: %i' % (tp,fp,tn,fn)
    
    print "Sensitivity: %f, Specificity: %f, Accuracy: %f" % (tp/(tp+fn), tn/(tn+fp), (tp+tn)/(tp+tn+fp+fn))
    
    #and add to total
    total_tp += tp
    total_tn += tn
    total_fp += fp
    total_fn += fn
    print total_tp, total_tn, total_fp, total_fn


# In[5]:

print "Sensitivity: %f, Specificity: %f, Accuracy: %f" % (total_tp/(total_tp+total_fn), 
                                                          total_tn/(total_tn+total_fp), 
                                                          (total_tp+total_tn)/(total_tp+total_tn+total_fp+total_fn))


# In[6]:

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

print model_pred
print true_labels

#plot ROC curves
fpr,tpr,_ = roc_curve(true_labels, scores)
auc = auc(fpr, tpr)

plt.plot(np.insert(fpr,0,0.), np.insert(tpr, 0,0.), label='Junction Counts (Area = %0.2f)' % auc, color="blue")
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel("False positive rate")
plt.ylabel("True positive rate")
plt.title("Deep Net Vorhersage")
plt.legend(loc="lower right")
#plt.show()
plt.savefig("ROC_jc.png", dpi=240)


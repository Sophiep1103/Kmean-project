import sys
import pandas as pd
import numpy as np
import math
import mykmeanssp

#extract the args from the command line to local variables.
def extract_args():
    if (len(sys.argv)!=4 and len(sys.argv)!=3):
        print("An Error Has Occurred")
        return -1
    if len(sys.argv)==3:
        k= -1 #no k given
        goal = sys.argv[1]
        file_name=sys.argv[2]
    if len(sys.argv)==4:
        k=sys.argv[1]
        goal = sys.argv[2]
        file_name=sys.argv[3]
    return k ,goal, file_name

#converts a file to list object
def fileToList(file_name):
    dataPoints = pd.read_csv(file_name, sep= ",", header= None)
    dataList= dataPoints.as_matrix()
    dataList=dataList.tolist()
    return dataList

#prints matrix
def printMatrix(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if j == len(matrix[0]) - 1:
                print(format(matrix[i][j], '.4f'), end="\n")
            else:
                print(format(matrix[i][j], '.4f'), end=",")



def main():
    input = extract_args()
    k=(int)(input[0])
    iter=300 #difault max iteration
    dataList = fileToList(input[2]) #a list of the data from the file given as an argument
    numPoints = len(dataList)
    pointLength= len(dataList[0])
    #checks which function should be execute (depands on the goal)
    if input[1]== "wam":
        try: #sending the datapoints to wam function, get weighted adj matrix
            W= mykmeanssp.wam(dataList, numPoints, pointLength)
            printMatrix(W)
        except:
            print("An Error Has Occurred")
            return 1
    elif input[1]=="ddg":
        try: #sending the datapoints to ddg function, get diagonal degree matrix
            D= mykmeanssp.ddg(dataList, numPoints, pointLength)
            printMatrix(D)
        except:
            print("An Error Has Occurred")
            return 1
    elif input[1]=="gl":
        try: #sending the datapoints to gl function, get graph Laplacian matrix
            L= mykmeanssp.gl(dataList, numPoints, pointLength)
            printMatrix(L)
        except:
            print("An Error Has Occurred")
            return 1
    elif input[1]=="jacobi":
        try: #sending the datapoints to jacobi function, get the eigenvalues and eigenvectors
            JACOBI= mykmeanssp.jacobi(dataList, numPoints, pointLength)
            printMatrix(JACOBI)
        except:
            print("An Error Has Occurred")
            return 1
    elif input[1]=="spk":
        try: #sending the datapoints to pre_spk function, get U back and then sending it to k_means++ algorithem
            U= mykmeanssp.pre_spk(dataList, k, numPoints, pointLength)
            #if k== -1 it should be calculated and and be the numbers of columns in U
            if k==-1:
                k= len(U[0])
            #calculate the solution to kmeans problem with U as points, starting with kmeans++ algorithm
            kmeans_pp_main(k, iter, U)
        except:
            print("An Error Has Occurred")
            return 1

#the main fuctions for kmeans_pp, the points are the rows in U
def kmeans_pp_main(k, iter, U):
    U_df = pd.DataFrame(U)
    indexes = [i for i in range(len(U))]
    U_df.insert(loc = 0, column='index', value= indexes) #add indexes to U
    epsilon=0.0
    means, means_indexes, numPoints, pointLenght = kmeans_pp(k, iter, epsilon, U, U_df)
    try:
        #sending the datapoints and initial centroids to spk function, get the final centroids
        kmeans_final_means = mykmeanssp.spk(k, iter, epsilon, U, numPoints, pointLenght, means)
        print_indexes(means_indexes)
        print_results(kmeans_final_means)
    except:
        print("An Error Has Occurred")
        return

#print the final means
def print_results(means):
    for point in means:
        for j in range(len(point)):
            if j == len(point) - 1:
                print(format(point[j], '.4f'), end="\n")
            else:
                print(format(point[j], '.4f'), end=",")

#print the means indexes
def print_indexes(indexes):
    for i in range(len(indexes)):
        if i == len(indexes) - 1:
            print(int(indexes[i]), end="\n")
        else:
            print(int(indexes[i]), end=",")

#calc the first means of the k_means algorithem
def kmeans_pp(k, iter, epsilon, U, U_df): 
    #number of points
    numPoints = len(U)
    #length of a point
    pointLenght = len(U[0])
    #array for the indexes of the points
    pointIndexes = [i for i in range(numPoints)]
    #get thr first random point index
    randomPointIndex = np.random.choice(pointIndexes)
    #initialize arrays for the means and their indexes
    means = [[0] for i in range(k)] 
    means_indexes = [0 for i in range(k)]
    #choose the first mean
    means[0] = U[randomPointIndex]
    means_indexes[0] = U_df.iloc[randomPointIndex, 0]
    num_means = 1
    #choosing the rest k-1 means to start with
    while (num_means < k):
        Dx = [0 for i in range(numPoints)] #create Dx
        prob = [0 for i in range(numPoints)] #create probs array
        calcDx(U, means, Dx, numPoints, num_means) #calculate Dx
        calcProbs(Dx, prob, numPoints) #calculate each probability for each point to be chosen
        #choose a point and update the relevent arrays
        chosen_element_index = np.random.choice(pointIndexes, p=prob) 
        means[num_means] = U[chosen_element_index]
        means_indexes[num_means] = U_df.iloc[chosen_element_index, 0]
        num_means += 1 #update the counter of the chosen means

    return means, means_indexes, numPoints, pointLenght

# calculate D(x) for all points
def calcDx(U, means, Dx, numPointsInner, num_means):
    for i in range(numPointsInner):
        dist = find_closest_mean(U[i], means[:num_means])
        Dx[i] = dist


#calcs the probailities array
def calcProbs(Dx, prob, numPointsInner):
    sum_dx = sum(Dx)
    for i in range(numPointsInner):
        prob[i] = Dx[i] / sum_dx  
    return prob

#calculates distance between point and mean
def distance(point, mean):
    dist = 0
    for i in range(len(point)):
        res = math.pow(point[i] - mean[i], 2)
        dist += res
    return math.sqrt(dist)


# calculate D(x) for a specific point, the minimal distance between point x and the means
def find_closest_mean(point, means):
    min_dist = sys.maxsize #initialize to the max number
    for i in range(len(means)):
        dist = distance(point, means[i])
        if dist < min_dist:
            min_dist = dist

    return min_dist


np.random.seed(0)
main()
    
 

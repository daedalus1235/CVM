import matplotlib.pyplot as plt
import csv


def main():
    files=['point', 'bond', 'square']

    plt.xlabel('Temperature')
    plt.ylabel('Composition of species 1')
    plt.title('Comparison of CVM Approximations for Ising Model')

    for filename in files:
        T=[]
        x=[]
        with open(filename+'.csv', 'r') as read:
            data = csv.reader(read, delimiter=',')
            for row in data:
                T.append(float(row[0]))
                x.append(float(row[1].strip(']').strip('[')))
        plt.plot(T,x, label=filename,alpha=0.5)

    plt.legend()
    plt.show()

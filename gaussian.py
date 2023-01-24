import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import statistics
import math

def plotGaussian(mMeans, mName, setN):

	x = np.arange(20, 130, 0.001)

	plt.plot(x, norm.pdf(x, mMeans[setN - 1][0][0], mMeans[setN - 1][0][1]), label=mName[setN - 1][0])
	plt.plot(x, norm.pdf(x, mMeans[setN - 1][1][0], mMeans[setN - 1][1][1]), label=mName[setN - 1][1])
	plt.plot(x, norm.pdf(x, mMeans[setN - 1][2][0], mMeans[setN - 1][2][1]), label=mName[setN - 1][2])

	plt.legend(title = 'Parameters')
	plt.xlabel('Means')
	plt.title("Set "+str(setN))
	plt.show()

def plotGaussianAll(mMeans, mName):

	x = np.arange(20, 130, 0.001)

	#Set 1
	plt.plot(x, norm.pdf(x, mMeans[0][0][0], mMeans[0][0][1]), label=mName[0][0])
	plt.plot(x, norm.pdf(x, mMeans[0][1][0], mMeans[0][1][1]), label=mName[0][1])
	plt.plot(x, norm.pdf(x, mMeans[0][2][0], mMeans[0][2][1]), label=mName[0][2])

	#Set 2
	plt.plot(x, norm.pdf(x, mMeans[1][0][0], mMeans[1][0][1]), label=mName[1][0])
	plt.plot(x, norm.pdf(x, mMeans[1][1][0], mMeans[1][1][1]), label=mName[1][1])
	plt.plot(x, norm.pdf(x, mMeans[1][2][0], mMeans[1][2][1]), label=mName[1][2])

	#Set 3
	plt.plot(x, norm.pdf(x, mMeans[2][0][0], mMeans[2][0][1]), label=mName[2][0])
	plt.plot(x, norm.pdf(x, mMeans[2][1][0], mMeans[2][1][1]), label=mName[2][1])
	plt.plot(x, norm.pdf(x, mMeans[2][2][0], mMeans[2][2][1]), label=mName[2][2])

	#Set 4
	plt.plot(x, norm.pdf(x, mMeans[3][0][0], mMeans[3][0][1]), label=mName[3][0])
	plt.plot(x, norm.pdf(x, mMeans[3][1][0], mMeans[3][1][1]), label=mName[3][1])
	plt.plot(x, norm.pdf(x, mMeans[3][2][0], mMeans[3][2][1]), label=mName[3][2])

	plt.legend(title='Parameters')
	plt.xlabel('Means')
	plt.title("All sets")
	plt.show()


def main():
	file = 'fingerprint/setScores.txt'
	f = open(file, "r")

	#Número de exemplos utilizados para testar
	nTeste = 9

	#Recebe número de sets e quantidade de variações das dimensões
	size = f.readline().split()
	n = int(size[0])
	m = int(size[1])

	#Inicializa matriz para dados coletados e nomes dos gráficos
	mData = np.empty((n, m, nTeste))
	mMeans = np.empty((n, m, 2))
	mName = [[0 for i in range(m)] for j in range(n)]

	#Número de sets
	for i in range(n):

		#Número de variações 
		for j in range(m):

			#Cabeçalhos dos dados - set e dimensão
			var = f.readline().split()
			setN = int(var[0])
			px = int(var[1])

			mName[i][j] = "Set "+str(setN)+" - "+str(px)+"px"

			#Gera matriz com dados
			for x in range(nTeste):
				mData[i][j][x] = float(f.readline().split()[0])

			mMeans[i][j][0] = statistics.mean(mData[i][j])
			mMeans[i][j][1] = statistics.stdev(mData[i][j])

			f.readline()

	print(mMeans)

	#Plot Graphics
	plotGaussian(mMeans, mName, 1)	#Set 1
	plotGaussian(mMeans, mName, 2)	#Set 2
	plotGaussian(mMeans, mName, 3)	#Set 3
	plotGaussian(mMeans, mName, 4)	#Set 4
	plotGaussianAll(mMeans, mName)	#Plot sll sets

	f.close()

if __name__ == "__main__":
    main()
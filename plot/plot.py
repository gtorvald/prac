import matplotlib.pyplot as plt

# data = []
# sumElem = 0
# for i in range(60):
# 	data.insert(i, float(input()))
# 	sumElem += data[i]

# print(sumElem / 60)

# plt.hist(data)
# plt.title("Распределение 1 - F, 26 кубитов, точность e = 0.001, 64 вычислительных узла, 60 запусков")
# plt.show()

qubits = [24, 25, 26, 27, 28]
data = [0.0044759, 0.0046625545, 0.004848487, 0.00503441, 0.0052206]
plt.bar(qubits, data)
plt.title("Среднее значение потерь точности")
plt.show()
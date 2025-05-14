# #Prime Number
import matplotlib.pyplot as plt

x = [1, 2, 3, 4]
y = [0.0595, 0.0316, 0.0224, 0.0336]
plt.figure(figsize=(10, 8))
plt.bar(x, y, color='skyblue', edgecolor='black')
plt.xlabel('Number of Processes')
plt.ylabel('Execution Time (seconds)')
plt.title('  Find Prime Number  ')
plt.xticks(x)
plt.tight_layout()
plt.show()


#Parallel Quick search 
import matplotlib.pyplot as plt
x2 = [1, 2, 3 ]
y2 = [0.000049, 0.000048, 0.0000400]
plt.figure(figsize=(10, 8))
plt.bar(x2, y2, color='skyblue', edgecolor='black')
plt.xlabel('Number of Processes')
plt.ylabel('Execution Time (seconds)')
plt.title('  Parallel Quick Search  ')
plt.xticks(x2)
plt.tight_layout()
plt.show()

#Radix Sort 
import matplotlib.pyplot as plt
x3 = [1, 2, 3 ,4]
y3 = [0.0399, 0.0219, 0.0285,0.0234]
plt.figure(figsize=(10, 8))
plt.bar(x3, y3, color='skyblue', edgecolor='black')
plt.xlabel('Number of Processes')
plt.ylabel('Execution Time (seconds)')
plt.title('  Radix Sort  ')
plt.xticks(x3)
plt.tight_layout()
plt.show()

#Bitonic Sort
import matplotlib.pyplot as plt
x4 = [1, 2, 4 ,8]
y4 = [0.03138, 0.04372, 0.08623,0.12016]
plt.figure(figsize=(10, 8))
plt.bar(x4, y4, color='skyblue', edgecolor='black')
plt.xlabel('Number of Processes')
plt.ylabel('Execution Time (seconds)')
plt.title('  Bitonic Sort ')
plt.xticks(x4)
plt.tight_layout()
plt.show()

#Sample Sort
import matplotlib.pyplot as plt
x4 = [1, 2, 3 ,4]
y4 = [0.0417, 0.0283, 0.0201,0.029]
plt.figure(figsize=(10, 8))
plt.bar(x4, y4, color='skyblue', edgecolor='black')
plt.xlabel('Number of Processes')
plt.ylabel('Execution Time (seconds)')
plt.title('  Sample Sort ')
plt.xticks(x4)
plt.tight_layout()
plt.show()


















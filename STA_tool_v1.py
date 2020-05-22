import math

# This table is copied from Probability and Statistical Inference
# written by Robert V. Hogg • Elliot A. Tanis • Dale L. Zimmerman.
NORMAL = [
    [0.5000, 0.5040, 0.5080, 0.5120, 0.5160, 0.5199, 0.5239, 0.5279, 0.5319, 0.5359],
    [0.5398, 0.5438, 0.5478, 0.5517, 0.5557, 0.5596, 0.5636, 0.5675, 0.5714, 0.5753],
    [0.5793, 0.5832, 0.5871, 0.5910, 0.5948, 0.5987, 0.6026, 0.6064, 0.6103, 0.6141],
    [0.6179, 0.6217, 0.6255, 0.6293, 0.6331, 0.6368, 0.6406, 0.6443, 0.6480, 0.6517],
    [0.6554, 0.6591, 0.6628, 0.6664, 0.6700, 0.6736, 0.6772, 0.6808, 0.6844, 0.6879],
    [0.6915, 0.6950, 0.6985, 0.7019, 0.7054, 0.7088, 0.7123, 0.7157, 0.7190, 0.7224],
    [0.7257, 0.7291, 0.7324, 0.7357, 0.7389, 0.7422, 0.7454, 0.7486, 0.7517, 0.7549],
    [0.7580, 0.7611, 0.7642, 0.7673, 0.7703, 0.7734, 0.7764, 0.7794, 0.7823, 0.7852],
    [0.7881, 0.7910, 0.7939, 0.7967, 0.7995, 0.8023, 0.8051, 0.8078, 0.8106, 0.8133],
    [0.8159, 0.8186, 0.8212, 0.8238, 0.8264, 0.8289, 0.8315, 0.8340, 0.8365, 0.8389],
    [0.8413, 0.8438, 0.8461, 0.8485, 0.8508, 0.8531, 0.8554, 0.8577, 0.8599, 0.8621],
    [0.8643, 0.8665, 0.8686, 0.8708, 0.8729, 0.8749, 0.8770, 0.8790, 0.8810, 0.8830],
    [0.8849, 0.8869, 0.8888, 0.8907, 0.8925, 0.8944, 0.8962, 0.8980, 0.8997, 0.9015],
    [0.9032, 0.9049, 0.9066, 0.9082, 0.9099, 0.9115, 0.9131, 0.9147, 0.9162, 0.9177],
    [0.9192, 0.9207, 0.9222, 0.9236, 0.9251, 0.9265, 0.9279, 0.9292, 0.9306, 0.9319],
    [0.9332, 0.9345, 0.9357, 0.9370, 0.9382, 0.9394, 0.9406, 0.9418, 0.9429, 0.9441],
    [0.9452, 0.9463, 0.9474, 0.9484, 0.9495, 0.9505, 0.9515, 0.9525, 0.9535, 0.9545],
    [0.9554, 0.9564, 0.9573, 0.9582, 0.9591, 0.9599, 0.9608, 0.9616, 0.9625, 0.9633],
    [0.9641, 0.9649, 0.9656, 0.9664, 0.9671, 0.9678, 0.9686, 0.9693, 0.9699, 0.9706],
    [0.9713, 0.9719, 0.9726, 0.9732, 0.9738, 0.9744, 0.9750, 0.9756, 0.9761, 0.9767],
    [0.9772, 0.9778, 0.9783, 0.9788, 0.9793, 0.9798, 0.9803, 0.9808, 0.9812, 0.9817],
    [0.9821, 0.9826, 0.9830, 0.9834, 0.9838, 0.9842, 0.9846, 0.9850, 0.9854, 0.9857],
    [0.9861, 0.9864, 0.9868, 0.9871, 0.9875, 0.9878, 0.9881, 0.9884, 0.9887, 0.9890],
    [0.9893, 0.9896, 0.9898, 0.9901, 0.9904, 0.9906, 0.9909, 0.9911, 0.9913, 0.9916],
    [0.9918, 0.9920, 0.9922, 0.9925, 0.9927, 0.9929, 0.9931, 0.9932, 0.9934, 0.9936],
    [0.9938, 0.9940, 0.9941, 0.9943, 0.9945, 0.9946, 0.9948, 0.9949, 0.9951, 0.9952],
    [0.9953, 0.9955, 0.9956, 0.9957, 0.9959, 0.9960, 0.9961, 0.9962, 0.9963, 0.9964],
    [0.9965, 0.9966, 0.9967, 0.9968, 0.9969, 0.9970, 0.9971, 0.9972, 0.9973, 0.9974],
    [0.9974, 0.9975, 0.9976, 0.9977, 0.9977, 0.9978, 0.9979, 0.9979, 0.9980, 0.9981],
    [0.9981, 0.9982, 0.9982, 0.9983, 0.9984, 0.9984, 0.9985, 0.9985, 0.9986, 0.9986],
    [0.9987, 0.9987, 0.9987, 0.9988, 0.9988, 0.9989, 0.9989, 0.9989, 0.9990, 0.9990]
    ]

def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n - 1)

# Calculate combinatorial number.
def combination(n, m):
    return factorial(n) / factorial(m) / factorial(n - m)

def binomial(n, p, lower_bound, upper_bound):
    answer = 0
    current = lower_bound
    while current <= upper_bound:
        answer += combination(n, current) * (p ** current) * ((1 - p) ** (n - current))
        current += 1
    return answer

def poisson(mean, lower_bound, upper_bound):
    answer = 0
    current = lower_bound
    while current <= upper_bound:
        answer += (mean ** current) * (math.e ** (mean * -1)) / factorial(current)
        current += 1
    return answer

# Change the general normal distribution into standard normal distribution.
def normal_transition(number, mean, variance):
    number = (number - mean) / math.sqrt(variance)
    number = round(number, 2) * 100
    if number > 3:
        number = 3
        print('Warning: This calculation is out of range. The outcome is approximated.')
    if number < -3:
        number = -3
        print('Warning: This calculation is out of range. The outcome is approximated.')
    row, column = divmod(number, 10)
    row, column = int(row), int(column)
    if number < 0:
        row *= -1
        return (1 - NORMAL[row][column])
    else:
        return NORMAL[row][column]

def normal(mean, variance, lower_bound, upper_bound):
    lower_bound = normal_transition(lower_bound, mean, variance)
    upper_bound = normal_transition(upper_bound, mean, variance)
    return upper_bound - lower_bound

# Prompt the user to enter certain variable and ensure the inputs are numbers.
def input_number(description, index):
    while True:
        number = input('Please enter {}: '.format(description)).replace(' ', '')
        if index == 0:
            try:
                number = int(number)
                return number
            except:
                print('Invalid input!')
        else:
            try:
                number = float(number)
                return number
            except:
                print('Invalid input!')

# 1 for binomial, 2 for Poisson, 3 for normal.
def choose(distribution):
    if distribution == '1':
        while True:
            n = input_number('n', 0)
            if n > 0:
                break
            else:
                print('n should be a positive integer.')
        while True:
            p = input_number('p', 1)
            if 0 < p < 1:
                break
            else:
                print('p should be a number between 0 and 1.')
        while True:
            lower_bound = input_number('lower bound', 0)
            if 0 <= lower_bound <= n:
                break
            else:
                print('lower_bound should be a positive integer less than n.')
        while True:
            upper_bound = input_number('upper bound', 0)
            if lower_bound <= upper_bound <= n:
                break
            else:
                print('upper bound should be an interger between lower bound and n.')
        answer = binomial(n, p, lower_bound, upper_bound)
        print('P({} <= X <= {}) = {}'.format(lower_bound, upper_bound, answer))
    elif distribution == '2':
        while True:
            mean = input_number('mean', 1)
            if mean > 0:
                break
            else:
                print('mean should be a positive number.')
        while True:
            lower_bound = input_number('lower bound', 0)
            if 0 <= lower_bound:
                break
            else:
                print('lower_bound should be a positive integer.')
        while True:
            print('If the upper bound is the positive infinite, please enter "-1"')
            upper_bound = input_number('upper bound', 0)
            if lower_bound <= upper_bound:
                break
            elif upper_bound == -1:
                break
            else:
                print('upper bound should be an interger no less than lower bound.')
        if upper_bound == -1:
            inverse_answer = poisson(mean, 0, lower_bound - 1)
            answer = 1 - inverse_answer
            print('P({} <= X) = {}'.format(lower_bound, answer))
        else:
            answer = poisson(mean, lower_bound, upper_bound)
            print('P({} <= X <= {}) = {}'.format(lower_bound, upper_bound, answer))
    elif distribution == '3':
        mean = input_number('mean', 1)
        while True:
            variance = input_number('variance', 1)
            if variance > 0:
                break
            else:
                print('variance should be a positive number.')
        lower_bound = input_number('lower bound', 1)
        while True:
            upper_bound = input_number('upper bound', 1)
            if lower_bound <= upper_bound:
                break
            else:
                print('upper bound should be no less than lower bound.')
        answer = normal(mean, variance, lower_bound, upper_bound)
        print('P({} <= X <= {}) = {}'.format(lower_bound, upper_bound, answer))

def main():
    # Some information.
    print('\n-----Welcome to STA tool 1.1!-----\n')
    print('This program is used for calculating the probability of \
binomial distribution, Poisson distribution, and normal distribution.\n')
    print('How to use this program:')
    print('First, enter the basic characteristics of the distribution \
X~b(n, p), X~Poisson(mean), and X~N(mean, variance).')
    print('Second, enter the lower bound and upper bound: \
P(lower bound <= X <= upper bound)')
    print('Then the program will tell you the probability.\nHope this will help you ^_^\n')
    go_on = True
    while go_on != '0':
        while True:
            print('----------------------------------------------------')
            distribution = input('Enter 1 to calculate with binomial distribution\n\
Enter 2 to calculate with Poisson distribution\nEnter 3 to calculate with normal distribution\n')
            if distribution == '1' or distribution == '2' or distribution == '3':
                break
            else:
                print('Invalid input.')
        choose(distribution)
        go_on = input('\nEnter 0 to quit the program, enter anything else to continue.\n')

main()
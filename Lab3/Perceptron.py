import numpy as np
import matplotlib.pyplot as plt
import os

def step_function(input_value):
    return 1 if input_value > 0 else 0

def initialize_weights(num_inputs):
    # Randomly initialize weights as numbers between -3 and 3
    # Add one weight for the bias
    return np.random.uniform(-3, 3, num_inputs + 1)

def initialize_expected(operation='OR'):
    # Initialize expected output for each example
    operation = operation.upper()
    if operation in ('UNKNOWN', 'OR'):
        return np.array([0, 1, 1, 1])
    if operation == 'AND':
        return np.array([0, 0, 0, 1])
    if operation == 'XOR':
        return np.array([0, 1, 1, 0])
    if operation == 'NAND':
        return np.array([1, 1, 1, 0])
    else:
        raise ValueError("Enter valid operation")

def initialize_input():
    # Initialize input examples. Each row is an example, each column is a feature
    # This is a 'batch' of 4 examples
    return np.array([[0, 0],
                     [1, 0],
                     [0, 1],
                     [1, 1]])

def perceptron(operation='OR'):
    INPUT = initialize_input()
    CORRECT_OUTPUT = initialize_expected(operation)

    NUM_EXAMPLES = INPUT.shape[0]
    NUM_INPUTS = INPUT.shape[1]

    W = initialize_weights(NUM_INPUTS)
    eta = 0.8  # Learning rate, should be between 0 and 1
    num_iterations = 0
    max_iter = 25 # Maximum number of iterations
    SaveOutput = []
    
    
    while num_iterations < max_iter:
        TOTAL_ERROR = np.zeros(NUM_EXAMPLES)
        
        for example in range(NUM_EXAMPLES):
            Total_signal = np.dot(INPUT[example], W[:NUM_INPUTS]) + W[NUM_INPUTS]  # Weighted sum + bias, input to the activation function
            Output = step_function(Total_signal) # Activation function
            Error = Output - CORRECT_OUTPUT[example] # Calcualte error between output and correct output
            TOTAL_ERROR[example] = abs(Error) # Store error for this example
            
            # Update weights according to the perceptron learning rule
            if Output > CORRECT_OUTPUT[example]:
                W[:NUM_INPUTS] -= abs(Error * INPUT[example] * eta)
                W[NUM_INPUTS] -= abs(Error * eta)  # Adjust bias
            elif Output < CORRECT_OUTPUT[example]:
                W[:NUM_INPUTS] += abs(Error * INPUT[example] * eta)
                W[NUM_INPUTS] += abs(Error * eta)  # Adjust bias
            
        num_iterations += 1
        SaveOutput.append(np.mean(TOTAL_ERROR))
        
        if np.sum(TOTAL_ERROR) == 0:
            break
    
    plt.plot(SaveOutput)
    plt.xlabel('Iterations')
    plt.ylabel('Mean Error')
    plt.title('Perceptron Training')
    plt.show()
    
    return W, num_iterations

# Run the perceptron function
operation = os.environ.get("OPERATION", "OR")
final_weights, iters = perceptron(operation=operation)
print("Operation:", operation.upper())
print("Final Weights:", final_weights)
print("Number of iterations:", iters)
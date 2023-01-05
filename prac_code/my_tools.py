import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

# few vertical plots
def few_v_plots(x_data, data, constr, fig_size=[25,30]):
    if not type(x_data) in [pd.Series, np.array, torch.Tensor, list]:
        raise Exception("type of x axis data can be one of these:"
                        "pd.Series, np.array, torch.Tensor, list."
                        "Not a {}".format(type(x_data)))
    
    x_data = pd.Series(x_data)
    
    if not x_data.dtype in [int, float]:
        raise Exception("elements of x axis data can have int or"
        "float type, not {}".format(x_data.dtype))
       
    size = x_data.size
    
    if type(data) == dict:
        data = pd.DataFrame(data)
    elif type(data) != pd.DataFrame:
        raise Exception("type of data can be pd.DataFrame or dictionary,"
        "not {}".format(type(data)))
    
    if type(constr) != list:
        raise Exception("type of construction can be list, not {}".format(type(constr)))
 
    num_axs = len(constr)
    fig, axs = plt.subplots(num_axs, 1, figsize = fig_size)
    
    if num_axs == 1:
        axs = [axs]
    
    for index in range(num_axs):    
        curve = constr[index]
        ax = axs[index]        
        if type(curve) == str:
            curve = [curve]
            
        if type(curve) == list:            
            for i in range(len(curve)):
                if type(curve[i]) != str:
                    raise Exceptrion("type of construction's elements"
                    "is string or list of strings, not {}".format(type(cirve[i])))
                elif not curve[i] in data.columns:
                    raise Exception("There is not column {} in data".format(curve[i]))
                else:    
                    ax.plot(x_data, data[curve[i]], label=curve[i])
                    ax.legend()
        else:
            raise Exception("type of construction's elements"
            "is string or list of strings, not {}".format(type(cirve[i])))
            
# plots few curves on one plot 
def plot_few_curves(x_datasets, y_datasets, curve_types, 
                    fig_size=[5,5], xlabel="", ylabel="", title=""):
    if not type(x_datasets) in [pd.Series, np.array, torch.Tensor, list]:
        raise Exception("type of x axis data can be one of these:"
                        "pd.Series, np.array, torch.Tensor, list."
                        "Not a {}".format(type(x_datasets)))
    
    if not type(y_datasets) in [pd.Series, np.array, torch.Tensor, list]:
        raise Exception("type of y axis data can be one of these:"
                        "pd.Series, np.array, torch.Tensor, list."
                        "Not a {}".format(type(y_datasets)))
    
    if type(x_datasets) != list:
        x_datasets = [x_datasets]
    
    if type(y_datasets) != list:
        y_datasets = [y_datasets]
    
    if type(curve_types) != list:
        raise Exception("curve types is a list, not {}".format(type(curve_types)))
    
    num_sets = len(x_datasets)
    if len(y_datasets) != num_sets:
        raise Exception("number of x datsets have to be equal number of y datasets")
    
    if len(curve_types) != num_sets:
        raise Exception("number of x datsets have to be equal number of curve types")
                        
    fig, ax = plt.subplots(1, 1, figsize=fig_size)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    for i in range(num_sets):
        if curve_types[i] == "plot":
            ax.plot(x_datasets[i], y_datasets[i])
        elif curve_types[i] == "scatter":
            ax.scatter(x_datasets[i], y_datasets[i])
        else:
            raise Exception("curve type have to be a plot or scatter,"
                            " not {}".format(curve_types[i]))
    
    return fig, ax


def optimize(f, der_f, 
             start_point, start_alpha, eps, max_step, 
             conj_method, step_method=None, c=None, cond="conj"):
    x = start_point
    conj_vec = -der_f(x)
    alpha = start_alpha
    step_count = 0
    history = [x] 
    
#     condition = либо conj (conjection) либо dist (distance)
#     это критерий останова по сопряженному вектору и по норме 
#     разности соседних шагов
    
    while step_count < max_step and \
    (cond!="conj" or np.linalg.norm(conj_vec) > eps) \
    and (cond!="dist" or step_count==0 or np.linalg.norm(x-prev_x) > eps):
#         новый альфа
        if step_method == "armijo":
            while f(x + alpha * conj_vec) > \
                f(x) + c * alpha * np.dot(conj_vec, der_f(x)):
                alpha /= 2

#         следующая точка
        next_x = x + alpha * conj_vec
        
#         новый сопряженный вектор
        y = der_f(next_x) - der_f(x)
        
        if conj_method == "DY":
            betta = np.dot(der_f(next_x),der_f(next_x))/np.dot(conj_vec, y)
        elif conj_method == "HZ":
            betta = np.dot(
            y-2*conj_vec*np.dot(y, y)/np.dot(conj_vec, y),
            der_f(next_x)/np.dot(conj_vec, y)
            )
        else:
            raise Exception("conj method can be DY or HZ, not {}".format(conj_method))
        conj_vec = -der_f(next_x) + betta * conj_vec
        
        step_count += 1
        
        prev_x = x
        x = next_x
        history.append(x)
#         следующий раз будем искать альфа начиная
#         с удвоенного предыдущего
        if step_method == "armijo":
            alpha = alpha*2
        
    return x, step_count, np.array(history)


def distance_curves(start_points, J, der_J, alpha, eps, max_iter, 
                    conj_method, step_method, c, fig_size, a, b, cond="conj"):    
    distances = []
    step_history = []
    
    num_lines = len(start_points)
    for start_point in start_points:
        min_x, step_number, history = optimize(J, der_J, 
                                  start_point, alpha, eps, max_iter, 
                                  conj_method, step_method, c, cond)
            
#         посмотрим как решение в ходе поиска по норме отличается от аналитического
        distance = (np.sum((history - (-a+b))**2, 1))**(1/2)
        distances.append(distance)
        step_history.append(range(step_number+1))
    
    plot_few_curves(step_history, distances, ["plot"]*num_lines, 
                    fig_size, xlabel=r'step', ylabel=r'$\|x_{step}-(-a+b)\|$',
                   title="Condition: {}, Conj. Method: {},"
                    " Step Method: {}, Max Step {}".format(cond, conj_method, step_method, max_iter))

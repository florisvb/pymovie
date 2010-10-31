import numpy as np

def compare(a, b=None, method=None):

    if type(a) is list:
        result = a[0]
        for i in range([1,len(a)]):
            if method == 'lighten':
                result = np.minimum(result,a[i])
            elif method == 'darken':
                result = np.maximum(result,a[i])
        return result
    elif b is not None:
        if method is 'lighten':
            result = np.minimum(a,b)
        elif method is 'darken':
            result = np.minimum(a,b)
            return lightened
    else:
        ValueError('please enter a valid a,b pair')
        
def darken(a, b=None):
    result = compare(a,b,method='darken')
def lighten(a, b=None):
    result = compare(a,b,method='lighten')
    
#def threshold()

def difference(a,b)

    diff = a-b
    diff<0









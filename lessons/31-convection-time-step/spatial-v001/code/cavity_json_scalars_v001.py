import numpy as np

def scalar_default(value):
 if isinstance(value,np.generic):return value.item()
 raise TypeError('Unsupported record value: '+type(value).__name__)

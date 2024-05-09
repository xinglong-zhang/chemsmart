import covalent as ct

local = ct.executor.LocalExecutor()

# Construct manageable tasks out of functions
# by adding the @covalent.electron decorator
@ct.electron(executor=local)
def add(x, y):
   return x + y

@ct.electron(executor=local)
def multiply(x, y):
   return x*y

@ct.electron(executor=local)
def divide(x, y):
   return x/y

# Construct the workflow by stitching together
# the electrons defined earlier in a function with
# the @covalent.lattice decorator
@ct.lattice
def workflow(x, y):
   r1 = add(x, y)
   r2 = [multiply(r1, y) for _ in range(4)]
   r3 = [divide(x, value) for value in r2]
   return r3

# Dispatch the workflow
dispatch_id = ct.dispatch(workflow)(1, 2)
result = ct.get_result(dispatch_id)
print(result)
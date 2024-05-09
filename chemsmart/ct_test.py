import covalent as ct

# Construct tasks as "electrons"
@ct.electron
def join_words(a, b):
    return ", ".join([a, b])

@ct.electron
def excitement(a):
    return f"{a}!"

# Construct a workflow of tasks
@ct.lattice
def simple_workflow(a, b):
    phrase = join_words(a, b)
    return excitement(phrase)

# Dispatch the workflow
dispatch_id = ct.dispatch(simple_workflow)("Hello", "World")


def raw_vector_sum(r):
    s = 0.
    for i in range(r.contents.size):
        s += r.contents.data[i]
    return s

def vector_sum(d):
    s = 0.
    for i in range(d.size()):
        s+=d[i]
    return s

def matrix_sum(m_input):
    s = 0
    for row in m_input:
        s += sum(row)
    return s

def tensor_sum(t_input):
    s = 0.
    for mx_input in t_input:
        s += matrix_sum(mx_input)
    return s

from libscientific.vector import * 
import random
from utils import (vector_sum,
                   raw_vector_sum)

def test_double_vector():
    # Create a random vector of 10 elements
    check_value_1 = 4.006685751243019
    random.seed(42)
    a = [random.random() for j in range(10)]
    d = DVector(a)
    
    # check the get function
    assert abs(d[1] - 0.0250) <= 1e-4

    # check for something wrong
    assert abs(vector_sum(d) - check_value_1) <= 1e-14

    # check the set function
    d[1] = -2
    # check for something wrong
    check_value_2 = 1.9816749960203515
    assert abs(vector_sum(d) - check_value_2) <= 1e-14
    
    # Convert to list
    dlst = d.tolist()
    # check for something wrong
    assert abs(sum(dlst) - check_value_2) <= 1e-14

    # Append value -123 at the end 
    dvector_append(d.dvect, -123)
    # check for something wrong
    assert abs(d[d.size()-1] - (-123)) <= 1e-3

    #Append value -123 at the end in a different way
    d.append(-123)
    # check for something wrong
    assert abs(d[d.size()-1] - (-123)) <= 1e-3

    # remove the value -2 at index 1
    dvector_remove_at(d.dvect, 1)
    # check for something wrong
    check_value_3 = -242.01832500397967
    assert abs(vector_sum(d) - check_value_3) <= 1e-14

    # Extend vector d with another vector
    b = [random.random() for j in range(4)]
    d.extend(b)
    # check for something wrong
    check_value_4 = -241.06895812070218
    assert abs(vector_sum(d) - check_value_4) <= 1e-14

    # Create a copy of d.d in q
    q = dvector_copy(d.dvect)
    # check for something wrong
    assert abs(raw_vector_sum(q)- check_value_4) <= 1e-14
    # delete q vector
    q_size = q.contents.size
    del_dvector(q)
    # check for something wrong
    assert q.contents.size != q_size
    
    # check if the double vector d have the value -123
    assert dvector_has_value(d.dvect, -123.0000) == 0

    # check if the double vector d have not the value -13
    assert dvector_has_value(d.dvect, -13.0000) == 1
    del d

def test_uivector():
    # Create a vector of integers from 0 to 6
    a = [0,1,2,3,4,5,6]
    u = UIVector(a)

    # check the get function
    assert abs(u[1] - 1) <= 1e-4

    # check for something worng
    assert abs(vector_sum(u) - sum(a)) <= 1e-14
    
    # set value 2 in poistion 1
    a[1] = 2
    u[1] = 2
    # check for something worng
    assert abs(vector_sum(u) - sum(a)) <= 1e-14

    # convert u to a list
    lst = u.tolist()
    # check for something worng
    assert abs(vector_sum(u) - sum(lst)) <= 1e-14

    # Add at the end the value 123
    a.append(123)
    uivector_append(u.uivect, 123)
    # check for something worng
    assert abs(vector_sum(u) - sum(a)) <= 1e-14

    # Append 123 at the end in a different way
    a.append(123)
    u.append(123)
    # check for something worng
    assert abs(vector_sum(u) - sum(a)) <= 1e-14

    # remove 2 at index 1")
    del a[1]
    uivector_remove_at(u.uivect, 1)
    # check for something worng
    assert abs(vector_sum(u) - sum(a)) <= 1e-14

    # Extend u with a vector b
    b = [7, 8, 9, 10]
    a.extend(b)
    u.extend(b)
    # check for something worng
    assert abs(vector_sum(u) - sum(a)) <= 1e-14

import random

from libscientific.vectlist import *


def test_vectlist():
    a = [[random.random() for j in range(10)] for i in range(3)]
    d = DVectorList(a)
    # d.debug()
    # print("get value")
    # vect.print_dvector(d[1])
    # print("set list")
    d[1] = [8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
    d.tolist()
    # print("print the list converted")
    # for item in dlst:
    #    print(item)

    # print("Append a list - high level")
    lst = [9, 8, 7, 6]
    d.append(lst)
    # print("Reappend ppend a list - low level")
    dvector_list_append(d.dvl, lst)
    # d.debug()

    # Test appending DVector object
    from libscientific.vector import DVector

    dv = DVector([1.0, 2.0, 3.0])
    d.append(dv)

    # Check len, iter, empty init, repr
    assert len(d) == 6
    vectors = list(d)
    assert len(vectors) == 6
    assert vectors[-1] == [1.0, 2.0, 3.0]

    d_empty = DVectorList()
    assert len(d_empty) == 0

    assert "DVectorList(" in repr(d)
    assert str(d) == repr(d)

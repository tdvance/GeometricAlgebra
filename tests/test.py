def asserteq(a, b=True):
    if a!=b:
        try:
            if abs(a - b) > 1E-10:
                print('a=', a)
                print('b=', b)
                print('diff=', abs(a-b))
                assert(a==b)
        except:
                print('a=', a)
                print('b=', b)
                assert(a==b)

            

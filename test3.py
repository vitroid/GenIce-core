from genice_core.topology import decompose_complex_path


if __name__ == "__main__":
    path = [1,2,3,4,5,6,7,8,9,2,10,11,12,13,14,15,12,16,1]
    print(path)
    for p in decompose_complex_path(path):
        print(p)

"""
  8 7 6
  9   5
1 2 3 4
  10 16
  11 12 13
     15 14

"""
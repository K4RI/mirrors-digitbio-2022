import sys
from CSCPattern import Pattern, PatternConfig, PartitionElement, pattern_parser
from CRPattern import Pattern as PatternCR
from diagram import init_diagram, add_intent, add_object, clean_flags

if __name__ == "__main__":
    theta = 1
    cfg = PatternConfig(theta)

    p1 = Pattern(instance='- - + +')
    print ('p1', p1)
    p2 = Pattern(instance='- - - +')
    print ('p2', p2)
    print ('p1 n p2', p1.intersect(p2))
    p3 = Pattern(instance='- - + +')
    print ('p3', p3)
    p4 = Pattern(instance='+ + - -')
    print ('p4', p4)

    print ('p1 == p2', p1 == p2)
    print ('p1 == p4', p1 == p4)

    pa = Pattern(dirty=False)
    pa.partition.add(PartitionElement(pattern_parser('- - -',1)))
    pa.partition.add(PartitionElement(pattern_parser('+ + -',4)))
    print ('pa', pa)
    pb = Pattern(dirty=False)
    pb.partition.add(PartitionElement(pattern_parser('- - -',1)))
    pb.partition.add(PartitionElement(pattern_parser('- - -',4)))
    print ('pb', pb)
    print ('pa == pb', pa==pb)

    p6 = PatternCR(instance='1 4 2 2', dirty=True)
    p7 = PatternCR(instance='3 3 2 2', dirty=True)
    print ('p6 =', p6)
    print ('p7 =', p7)
    print ('p6 n p7 =', p6.intersect(p7))
    print ('p6 == p7', p6 == p7)
    print ('p6 <= p7', p6 <= p7)

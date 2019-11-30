
class Node(object):

    def __init__(self,id: int, x: float, y: float, value: float = 0,BC: bool = False):
        self.__id = id
        self.__x = x
        self.__y = y
        self.__value = value
        self.__BC = BC

    @property
    def id(self):
        return self.__id

    @id.setter
    def id(self, new_value: float):
        self.__id = new_value

    @property
    def x(self):
        return self.__x

    @x.setter
    def x(self, new_value: float):
        self.__x = new_value

    @property
    def y(self):
        return self.__y

    @y.setter
    def y(self, new_value: float):
        self.__y = new_value

    @property
    def value(self):
        return self.__y

    @value.setter
    def value(self, new_value: float):
        self.__y = new_value

    @property
    def BC(self):
        return self.__BC

    @BC.setter
    def BC(self, new_value: float):
        self.__BC = new_value

    def __str__(self):
        return "N(id:{}, x:{}, y:{}, value:{}, BC:{})".format(self.__id, self.__x, self.__y, self.__value, self.__BC)

    def __repr__(self):
        return str(self)


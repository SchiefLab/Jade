

class Tester:
	def __init__(self):
		pass
	
	def my_method(self, *args, **kwargs):
		print repr(args)
		print repr(kwargs)


class Tester2(Tester):
	def __init__(self):
		pass
	
	def my_method(self, v):
		print repr(v)

if __name__ == "__main__":
	test1 = Tester()
	test1.my_method()
	print "\n\n"

	test2 = Tester2()
	test2.my_method("So, did it work?")


	print("Yes indeedy!")	

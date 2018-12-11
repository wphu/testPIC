default:

clean:
	make -C timer							clean
	make -C timer-class-variables			clean
	make -C ref								clean
	make -C ref--function					clean
	make -C ref-class						clean
	make -C ref-class-virtual-function		clean
	make -C ref-class-aos					clean
	make -C ref-class-soa					clean

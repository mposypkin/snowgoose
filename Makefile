all dep clean indent tests::
	cd libjson && make $@ && cd .. \\
	cd common && make $@ && cd .. \\
	cd interval && make $@ && cd

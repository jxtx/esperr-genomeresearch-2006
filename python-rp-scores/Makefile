inplace: setup
	cp target/rp/standard_model.so rp

setup:
	mkdir -p target
	python setup.py install --install-lib=target --install-scripts=target

test:
	( cd target; python test_standard_model.py )

clean:
	python setup.py clean

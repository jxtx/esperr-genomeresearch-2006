inplace: setup
	cp target/rp/models/standard.so rp/models
	cp target/rp/mapping_helper.so rp

setup:
	mkdir -p target
	python setup.py install --install-lib=target --install-scripts=target

test:
	( cd target; python test_standard_model.py )

clean:
	python setup.py clean

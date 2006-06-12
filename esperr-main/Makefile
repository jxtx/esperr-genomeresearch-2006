# Just wraps setup.py

inplace: setup
	cp target/rp/models/standard.so rp/models
	#cp target/rp/models/simple_periodic.so rp/models
	#cp target/rp/models/complex_periodic.so rp/models
	#cp target/rp/models/tree.so rp/models
	#cp target/rp/models/tree_pruned_1.so rp/models
	cp target/rp/models/tree_pruned.so rp/models
	#cp target/rp/models/gill_pst.so rp/models
	#cp target/rp/mapping_helper.so rp

setup:
	mkdir -p target
	python setup.py build --debug install --install-lib=target --install-scripts=target

test: inplace
	python test_standard_model.py

clean:
	python setup.py clean
	rm -rf build

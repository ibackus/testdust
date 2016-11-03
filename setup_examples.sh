#!/bin/bash
# Check for overwriting
if [ -d "examples/" ]; then
    while true; do
        read -p "Overwrite your copies of examples [y/n]? " yn
        case $yn in
            [Yy]* ) break;;
            [Nn]* ) echo "skipping"; exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done
fi
# Copy examples
rsync -avI src/examples/ examples/
#rm examples/README.md
echo "Success.  See examples in examples/ folder"


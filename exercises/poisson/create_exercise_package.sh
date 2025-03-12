set -e

PACKAGE_FOLDER_NAME=poisson_exercise

rm -rf $PACKAGE_FOLDER_NAME
mkdir -p $PACKAGE_FOLDER_NAME/C++
mkdir -p $PACKAGE_FOLDER_NAME/Fortran

cp -r C++/cpu_original $PACKAGE_FOLDER_NAME/C++
cp -r Fortran/cpu_original $PACKAGE_FOLDER_NAME/Fortran

cp requirements.txt $PACKAGE_FOLDER_NAME
cp env-archer2.sh $PACKAGE_FOLDER_NAME
cp check_output.py $PACKAGE_FOLDER_NAME


tar -zcvf $PACKAGE_FOLDER_NAME.tar.gz $PACKAGE_FOLDER_NAME
rm -r $PACKAGE_FOLDER_NAME

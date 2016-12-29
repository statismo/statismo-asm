#include <itkDirectory.h>
#include <itkStandardMeshRepresenter.h>
#include <ASM/itkASMNormalDirectionFeatureExtractor.h>
#include <ASM/itkASMGaussianGradientImagePreprocessor.h>
#include <ASM/itkActiveShapeModelBuilder.h>
#include <ASM/itkASMRandomPointIdSampler.h>
#include <ASM/itkDeferredItems.h>

typedef itk::Mesh<float, 3> MeshType;
typedef itk::Image<float, 3> ImageType;
typedef itk::ActiveShapeModel<MeshType, ImageType> ActiveShapeModelType;
typedef itk::ActiveShapeModelBuilder<MeshType, ImageType> ASMBuilderType;
typedef itk::MeshFileReader<MeshType> MeshReaderType;

typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
typedef std::vector<std::string> StringVectorType;


void loadDirectory (std::string dir, std::vector<std::string> &files, const std::string& extension=".*") {
    itk::Directory::Pointer directory = itk::Directory::New();
    directory->Load(dir.c_str());

    for (unsigned i = 0; i < directory->GetNumberOfFiles(); i++) {
        const char* filename = directory->GetFile(i);
        std::string sfilename(filename);
        if (extension == ".*" || sfilename.rfind(extension) == sfilename.length() - extension.length())
            files.push_back(filename);
    }

}

std::string basename(std::string fullname) {
    size_t lastindex = fullname.find_last_of(".");
    return fullname.substr(0, lastindex);
}


int main(int argc, char *argv[]) {


    std::string imagesDirectory("/tmp/ulna/i");
    std::string meshesDirectory("/tmp/ulna/m");

    StringVectorType imageFiles;
    loadDirectory(imagesDirectory, imageFiles, ".nii");

    StringVectorType meshFiles;
    loadDirectory(meshesDirectory, meshFiles, ".vtk");


    // we need pairs of images and meshes. The matching is done by file base name (i.e., without the extension)
    ASMBuilderType::DataPairsType data;

    std::map<std::string, std::string> imageNames;

    for (StringVectorType::const_iterator imageFile=imageFiles.begin(); imageFile != imageFiles.end(); imageFile++) {
        std::string fullpath = (std::string(imagesDirectory) + "/") + *imageFile;
        imageNames[basename(*imageFile)] = fullpath;
    }

    // matchmaking
    for (StringVectorType::const_iterator meshFile=meshFiles.begin(); meshFile != meshFiles.end(); meshFile++) {
        std::string key = basename(*meshFile);
        std::map<std::string, std::string>::iterator it = imageNames.find(key);
        if (it != imageNames.end()) {
            std::cout << "Using: " << key << std::endl;
            std::string meshPath = (std::string(meshesDirectory) + "/") + *meshFile;
            std::string imagePath = it->second;

            itk::DeferredImage<ImageType>* image = new itk::DeferredImage<ImageType>(imagePath);
            itk::DeferredMesh<MeshType>* mesh = new itk::DeferredMesh<MeshType>(meshPath);

            ASMBuilderType::DataPairType pair(key, std::make_pair(mesh, image));
            data.push_back(pair);
        }
    }

    MeshReaderType::Pointer reader = MeshReaderType::New();
    reader->SetFileName("/tmp/ulna/m/vsd-7.vtk");
    reader->Update();
    MeshType::Pointer refMesh = reader->GetOutput();

    RepresenterType::Pointer representer = RepresenterType::New();
    representer->SetReference(refMesh);

    // sampler to determine where to create profiles
    itk::ASMRandomPointIdSampler<MeshType>::Pointer pointIdSampler = itk::ASMRandomPointIdSampler<MeshType>::New();
    pointIdSampler->SetMesh(refMesh);
    pointIdSampler->SetNumberOfPoints(2600);
    pointIdSampler->Update();

    // image preprocessor
    itk::ASMGaussianGradientImagePreprocessor<MeshType, ImageType>* preprocessor = new itk::ASMGaussianGradientImagePreprocessor<MeshType, ImageType>(2.0);

    // feature extractor
    itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::Pointer profileSampler = itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::New();
    profileSampler->SetNumberOfPoints(7);
    profileSampler->SetPointSpacing(1);
    itk::ASMNormalDirectionFeatureExtractor<MeshType, ImageType>* featureExtractor = new itk::ASMNormalDirectionFeatureExtractor<MeshType, ImageType>(profileSampler, 1.0);

    // build!
    ASMBuilderType::Pointer asmBuilder = ASMBuilderType::New();
    ActiveShapeModelType::Pointer newModel = asmBuilder->BuildNewModel(representer, data, featureExtractor, preprocessor, pointIdSampler);
    newModel->Save("/tmp/asm-c.h5");

    std::cout << "DONE" <<std::endl;

    return 0;
}


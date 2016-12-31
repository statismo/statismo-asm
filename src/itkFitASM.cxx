#include <iostream>
#include <itkMesh.h>
#include <itkImage.h>
#include <itkImageFileReader.h>

#include <itkStandardMeshRepresenter.h>
#include "ASM/ASMFitting.h"
#include <itkEuler3DTransform.h>
#include "ASM/itkASMNormalDirectionPointSampler.h"
#include "ASM/itkASMNormalDirectionFeatureExtractor.h"
#include "ASM/itkASMGaussianGradientImagePreprocessor.h"
#include "ASM/itkASMIdentityImagePreprocessor.h"
#include "ASM/itkActiveShapeModel.h"
#include "ASM/itkASMFitting.h"
#include "itkTimeProbe.h"
#include "ASM/itkASMOps.h"

typedef itk::Mesh<float, 3> MeshType;
typedef itk::Image<float, 3> ImageType;
typedef itk::ActiveShapeModel<MeshType, ImageType> ActiveShapeModelType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;

typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
typedef itk::ASMMeshOperations::RigidTransformType RigidTransformType;
typedef itk::ASMMeshOperations::RigidTransformPointerType RigidTransformPointerType;

typedef RepresenterType::PointType PointType;
typedef itk::StatisticalModel<MeshType> StatisticalModelType;
typedef itk::Euler3DTransform< float > TransformType;

typedef itk::ASMFitting<MeshType, ImageType> FittingType;
typedef itk::ASMFittingStep<MeshType, ImageType> FittingStepType;
typedef itk::ASMFittingResult<MeshType, ImageType> FittingResultType;


// FIXME: these conversions have to go.
typedef vnl_vector<statismo::ScalarType> VnlVectorType;
VnlVectorType toVnlVector(const statismo::VectorType& v) {
    return VnlVectorType(v.data(), v.rows());
}
statismo::VectorType fromVnlVector(const VnlVectorType& v) {
    return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());

}

int main(int argc, char *argv[]) {


    if (argc < 4) {
        std::cerr << "usage: " << argv[0] << " modelname targetImage outputDir" << std::endl;
        exit(-1);
    }

    std::string modelname(argv[1]);
    std::string targetImageName(argv[2]);
    std::string outputDir(argv[3]);


    std::cout << "Initializing..." << std::endl;
    // FIXME: these should go somewhere less "intrusive"
    statismo::ASMFeatureExtractorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMNormalDirectionFeatureExtractorFactory<MeshType, ImageType>::GetInstance());
    statismo::ASMImagePreprocessorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMGaussianGradientImagePreprocessorFactory<MeshType, ImageType>::GetInstance());
    statismo::ASMImagePreprocessorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMIdentityImagePreprocessorFactory<MeshType, ImageType>::GetInstance());

    RepresenterType::Pointer representer = RepresenterType::New();


    ActiveShapeModelType::Pointer aModel = ActiveShapeModelType::New();
    aModel->Load(representer,  modelname.c_str());

    itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::Pointer fitSampler = itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::New();
    fitSampler->SetNumberOfPoints(25);
    fitSampler->SetPointSpacing(1);

    statismo::ASMFittingConfiguration fitConfig(3,5,3);

    RigidTransformType::Pointer currentTransform = RigidTransformType::New();
    currentTransform->SetIdentity();




    ImageReaderType::Pointer imgReader = ImageReaderType::New();
    imgReader->SetFileName(targetImageName.c_str());

    imgReader->Update();
    ImageType::Pointer image = imgReader->GetOutput();
    statismo::ASMPreprocessedImage<MeshType> *pimage = aModel->GetstatismoImplObj()->GetImagePreprocessor()->Preprocess(image);

    // just for testing
    aModel->SetStatisticalModel(aModel->GetStatisticalModel());

    statismo::VectorType coeffs = statismo::VectorType::Zero(aModel->GetStatisticalModel()->GetNumberOfPrincipalComponents());

    FittingStepType::Pointer fittingStep = FittingStepType::New();
    fittingStep->SetModel(aModel);
    fittingStep->SetTarget(pimage);
    fittingStep->SetSampler(FittingStepType::SamplerPointerType(fitSampler.GetPointer()));
    fittingStep->SetConfiguration(fitConfig);

    std::cout << "Initialization done." << std::endl;

    for (int i =1; i <= 10; ++i) {
        itk::TimeProbe clock;
        clock.Start();
        std::cout << "iteration: " << i << std::endl;
        fittingStep->SetCoefficients(coeffs);
        fittingStep->SetRigidTransformation(currentTransform);
        fittingStep->Update();
        FittingResultType::Pointer result = fittingStep->GetOutput();
        if (!result->IsValid()) {
            std::cout << "invalid result, aborting " <<std::endl;
            exit(42);
        }
        coeffs = fromVnlVector(result->GetCoefficients());
        currentTransform = result->GetRigidTransformation();
        std::cout << "coeffs (adj)" << toVnlVector(coeffs) << std::endl;
        clock.Stop();

        double elapsed = clock.GetMean();
        std::cout << "Elapsed " << elapsed << std::endl;
        if (currentTransform) {
            std::cout << "Writing result of iteration " << i << std::endl;
            itk::MeshFileWriter<MeshType>::Pointer writer = itk::MeshFileWriter<MeshType>::New();
            std::stringstream filename;
            filename << outputDir << "asmfit-iteration-" << i << ".vtk";
            writer->SetFileName(filename.str());
            MeshType::Pointer ms = result->GetMesh();
            writer->SetInput(ms);
            writer->Update();
        }
    }
    return 0;
}


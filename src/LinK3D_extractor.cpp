#include "LinK3D_extractor.h"

using namespace std;

namespace LinK3D_SLAM
{
    LinK3D_Extractor::LinK3D_Extractor(int nScans, float minimumRange, float XYDisTh, 
                        float xDisTh, float yDisTh, int ptNumTh, int scanNumTh, int scoreTh):
                        mNScans(nScans), mMinimumRange(minimumRange), mXYDisTh(XYDisTh),
                        mXDisTh(xDisTh), mYDisTh(yDisTh), mPtNumTh(ptNumTh), mScanNumTh(scanNumTh), mScoreTh(scoreTh){}

    /// @brief 剔除离载体较近的点云
    /// @param cloudIn 
    /// @param cloudOut 
    void LinK3D_Extractor::removeClosedPoint(pcl::PointCloud<pcl::PointXYZ> &cloudIn, pcl::PointCloud<pcl::PointXYZ> &cloudOut)
    {
        std::vector<int> indices;
        pcl::removeNaNFromPointCloud(cloudIn, cloudIn, indices); // PCL库函数 完成非法值点云的提取
        // 根据剔除非法值后的点云数量重新分配内存
        if(&cloudIn != &cloudOut){
            cloudOut.header = cloudIn.header;
            cloudOut.points.resize(cloudIn.points.size());
        }

        size_t j = 0;

        // mMinimumRange 可以理解为一个人为设定的“盲区”，离载体原点太近就剔除该点云
        for(size_t i = 0; i < cloudIn.points.size(); ++i)
        {
            if(cloudIn.points[i].x * cloudIn.points[i].x + cloudIn.points[i].y * cloudIn.points[i].y + cloudIn.points[i].z * cloudIn.points[i].z < mMinimumRange * mMinimumRange)
                continue;
            cloudOut.points[j] = cloudIn.points[i];
            j++;
        }

        if(j != cloudIn.points.size())
        {
            cloudOut.points.resize(j);
        }

        cloudOut.height = 1;
        cloudOut.width = static_cast<uint32_t>(j);
        cloudOut.is_dense = true;
    }

    /// @brief 边缘点提取
    /// Step 1 移除距离载体太近的点云
    /// Step 2 计算点云基础信息——包括scanID、水平角等
    /// Step 3：计算点云平滑度(或者说曲率)
    /// @param laserCloudIn 输入原始点云
    /// @param scanCloud 输出带有平滑度、scanID等信息的点云
    void LinK3D_Extractor::getEdgePoint(pcl::PointCloud<pcl::PointXYZ> &laserCloudIn, MatPt &scanCloud)
    {
        // 根据scan数量分配内存
        std::vector<int> scanStartInd(mNScans, 0);
        std::vector<int> scanEndInd(mNScans, 0);       

        // Step 1 移除距离载体太近的点云
        removeClosedPoint(laserCloudIn, laserCloudIn);

        // Step 2 计算点云基础信息——包括scanID、水平角等
        // 这一部分和A-LOAM一样，可参考其相关注释
        int cloudSize = laserCloudIn.points.size();
        
        // 计算起点和终点的角度，取反从原来的顺时针旋转变为逆时针旋转
        // 首先atan2的范围是[-pi, pi] 为了保证起点和终点的角度相差大致为2pi，所以这里加上2pi是为了和实际保持一致
        float startOri = -atan2(laserCloudIn.points[0].y, laserCloudIn.points[0].x);
        float endOri = -atan2(laserCloudIn.points[cloudSize - 1].y,
                          laserCloudIn.points[cloudSize - 1].x) + 2 * M_PI;
        
        // 理想情况是起点和终点的角度相差大致为2pi，这些对一些特殊情况进行补偿
        if (endOri - startOri > 3 * M_PI)
        {
            endOri -= 2 * M_PI;
        }
        else if (endOri - startOri < M_PI)
        {
            endOri += 2 * M_PI;
        }
        
        // halfPassed表示点云的水平角是否超过了pi
        bool halfPassed = false;
        int count = cloudSize;
        pcl::PointXYZI point;
        std::vector<pcl::PointCloud<pcl::PointXYZI>> laserCloudScans(mNScans);

        for (int i = 0; i < cloudSize; i++)
        {
            point.x = laserCloudIn.points[i].x;
            point.y = laserCloudIn.points[i].y;
            point.z = laserCloudIn.points[i].z;

            // 计算点云的俯仰角，进而确定点云在那一条scan上
            float angle = atan(point.z / sqrt(point.x * point.x + point.y * point.y)) * 180 / M_PI;
            int scanID = 0;

            // scan线通常不是均匀分布的
            if (mNScans == 16)
            {
                scanID = int((angle + 15) / 2 + 0.5);
                if (scanID > (mNScans - 1) || scanID < 0)
                {
                    count--;
                    continue;
                }
            }
            else if (mNScans == 32)
            {
                scanID = int((angle + 92.0/3.0) * 3.0 / 4.0);
                if (scanID > (mNScans - 1) || scanID < 0)
                {
                    count--;
                    continue;
                }
            }
            else if (mNScans == 64)
            {   
                if (angle >= -8.83)
                    scanID = int((2 - angle) * 3.0 + 0.5);
                else
                    scanID = mNScans / 2 + int((-8.83 - angle) * 2.0 + 0.5);

                if (angle > 2 || angle < -24.33 || scanID > 50 || scanID < 0)
                {
                    count--;
                    continue;
                }
            }
            else
            {
                printf("wrong scan number\n");
            }
            
            // 计算点云的水平角，可以大致判断点云在空间中的位置 并且根据和起始角的关系得到当前点云的时间戳
            float ori = -atan2(point.y, point.x);
            if (!halfPassed)
            { 
                if (ori < startOri - M_PI / 2)
                {
                    ori += 2 * M_PI;
                }
                else if (ori > startOri + M_PI * 3 / 2)
                {
                    ori -= 2 * M_PI;
                }

                if (ori - startOri > M_PI)
                {
                    halfPassed = true;
                }
            }
            else
            {
                ori += 2 * M_PI;
                if (ori < endOri - M_PI * 3 / 2)
                {
                    ori += 2 * M_PI;
                }
                else if (ori > endOri + M_PI / 2)
                {
                    ori -= 2 * M_PI;
                }
            }

            // 这里的话 和A-LOAM的处理有些许区别，直接保存了点云的水平角度，并没有直接计算点云时间戳
            // TODO：
            point.intensity = ori;
            // 按照scanID保存点云
            laserCloudScans[scanID].points.push_back(point);
        }

        size_t scanSize = laserCloudScans.size();
        scanCloud.resize(scanSize);
        cloudSize = count;
        
        // Step 3：计算点云平滑度(或者说曲率)
        for(int i = 0; i < mNScans; i++)
        {
            int laserCloudScansSize = laserCloudScans[i].size();
            if(laserCloudScansSize >= 15)
            {
                // 根据当前点云的前后五个点
                for(int j = 5; j < laserCloudScansSize - 5; j++)
                {
                    float diffX = laserCloudScans[i].points[j - 5].x + laserCloudScans[i].points[j - 4].x
                                + laserCloudScans[i].points[j - 3].x + laserCloudScans[i].points[j - 2].x
                                + laserCloudScans[i].points[j - 1].x - 10 * laserCloudScans[i].points[j].x
                                + laserCloudScans[i].points[j + 1].x + laserCloudScans[i].points[j + 2].x
                                + laserCloudScans[i].points[j + 3].x + laserCloudScans[i].points[j + 4].x
                                + laserCloudScans[i].points[j + 5].x;
                    float diffY = laserCloudScans[i].points[j - 5].y + laserCloudScans[i].points[j - 4].y
                                + laserCloudScans[i].points[j - 3].y + laserCloudScans[i].points[j - 2].y
                                + laserCloudScans[i].points[j - 1].y - 10 * laserCloudScans[i].points[j].y
                                + laserCloudScans[i].points[j + 1].y + laserCloudScans[i].points[j + 2].y
                                + laserCloudScans[i].points[j + 3].y + laserCloudScans[i].points[j + 4].y
                                + laserCloudScans[i].points[j + 5].y;
                    float diffZ = laserCloudScans[i].points[j - 5].z + laserCloudScans[i].points[j - 4].z
                                + laserCloudScans[i].points[j - 3].z + laserCloudScans[i].points[j - 2].z
                                + laserCloudScans[i].points[j - 1].z - 10 * laserCloudScans[i].points[j].z
                                + laserCloudScans[i].points[j + 1].z + laserCloudScans[i].points[j + 2].z
                                + laserCloudScans[i].points[j + 3].z + laserCloudScans[i].points[j + 4].z
                                + laserCloudScans[i].points[j + 5].z;

                    float smoothness = diffX * diffX + diffY * diffY + diffZ * diffZ;
                    if(smoothness > 10 && smoothness < 20000)  
                    {
                        float ori = laserCloudScans[i].points[j].intensity;
                        
                        // ！！！
                        PointXYZSCA tmpPt;
                        tmpPt.x = laserCloudScans[i].points[j].x;
                        tmpPt.y = laserCloudScans[i].points[j].y;
                        tmpPt.z = laserCloudScans[i].points[j].z;
                        tmpPt.scan = i;
                        tmpPt.smoothness = smoothness;
                        tmpPt.angel = ori; 
                        scanCloud[i].emplace_back(tmpPt);
                    }
                }
            }
        }
       
    }

    /// @brief 根据每个点云的水平角确定其扇形区域位置
    /// @param scanCloud 输入：根据scanID存储的点云
    /// @param areaCloud 输出根据扇形区域存储点云
    void LinK3D_Extractor::divideArea(MatPt &scanCloud, MatPt &areaCloud)
    {
        int scanSize = scanCloud.size();
        if(scanSize == 0){
            return;
        }

        // 扇形数量 Nsect = 120
        areaCloud.resize(120);
        for(int i = 0; i < scanSize; i++) 
        {
            int pointSize = scanCloud[i].size();
            for(int j = 0; j < pointSize; j++)
            {                
                int areaNum = 0;
                float angel = scanCloud[i][j].angel;
                
                // 根据水平角 确定 在哪个扇形区域里
                if(angel > 0 && angel < 2 * M_PI)
                {
                    areaNum = std::floor((angel / (2 * M_PI)) * 120);
                }   
                else if(angel > 2 * M_PI)
                {
                    areaNum = std::floor(((angel - 2 * M_PI) / (2 * M_PI)) * 120);
                }
                else if(angel < 0)
                {
                    areaNum = std::floor(((angel + 2 * M_PI) / (2 * M_PI)) * 120);
                }
                areaCloud[areaNum].push_back(scanCloud[i][j]);
            }
        }
    }

    /// @brief 计算当前聚类中在XoY平面到原点的均值距离
    /// @param oneCluster 
    /// @return 
    float LinK3D_Extractor::computeClusterMeans(std::vector<PointXYZSCA> &oneCluster)
    {
        if(oneCluster.empty()){
            return 0;
        }

        float distSum = 0;
        int ptCnt = oneCluster.size();

        for(int i = 0; i < ptCnt; i++)
        {
            distSum += distXY(oneCluster[i]);
        }

        return (distSum/ptCnt);
    }

    /// @brief 分别计算当前聚类中所有点在x y方向上的位置均值
    /// @param oneCluster 
    /// @param xyMeans 
    void LinK3D_Extractor::computeXYMeans(std::vector<PointXYZSCA> &oneCluster, std::pair<float, float> &xyMeans)
    {
        if(oneCluster.empty()){
            return;
        }

        int ptCnt = oneCluster.size();
        float xSum = 0;
        float ySum = 0;

        for(int ptNum = 0; ptNum < ptCnt; ptNum++)
        {
            xSum += oneCluster[ptNum].x;
            ySum += oneCluster[ptNum].y;
        }

        float xMeans = xSum/ptCnt;
        float yMeans = ySum/ptCnt;
        xyMeans = std::make_pair(xMeans, yMeans);
    }

    /// @brief 点云聚类
    /// Step 1 所有扇形区域分别完成点云聚类
    /// Step 2 找到哪些簇需要合并
    /// Step 3 合并聚类簇
    /// @param areaCloud 输入按照扇形区域划分的点云
    /// @param clustered 输出聚类后的点云簇
    void LinK3D_Extractor::computeCluster(const MatPt &areaCloud, MatPt &clustered)
    {    
        int areaSize = areaCloud.size();
        if(areaSize == 0){
            return;
        }

        MatPt tmpClustered;
        PointXYZSCA curvPt;
        std::vector<PointXYZSCA> dummy(1, curvPt); 

        //Cluster for each sector area
        // Step 1 所有扇形区域分别完成点云聚类
        // 每个扇形区域分别聚类
        // 遍历每个扇形区域
        for(int i = 0; i < areaSize; i++)
        {
            if(areaCloud[i].size() < 6)  
                continue;

            int ptCnt = areaCloud[i].size();        
            MatPt curAreaCluster(1, dummy);
            // areaCloud[i][0]当前扇形区域中的第一个点，用于作为当前区域中第一个聚类点中的第一个点
            curAreaCluster[0][0] = areaCloud[i][0];

            // 遍历当前扇形区域中的所有点云
            for(int j = 1; j < ptCnt; j++)
            {
                // 当前扇形区域中的聚类数量
                int clusterCnt = curAreaCluster.size();

                // 遍历所有聚类
                for(int k = 0; k < clusterCnt; k++)
                {
                    // 计算当前聚类中在XoY平面到原点的均值距离
                    float means = computeClusterMeans(curAreaCluster[k]);
                    std::pair<float, float> xyMeans;
                    // 分别计算当前聚类中所有点在x y方向上的位置均值
                    computeXYMeans(curAreaCluster[k], xyMeans);
                    
                    PointXYZSCA tmpPt = areaCloud[i][j];
                    
                    // 在“均值点”范围内的点被聚类在一起
                    if(abs(distXY(tmpPt) - means) < mXYDisTh && abs(xyMeans.first - tmpPt.x) < mXDisTh && abs(xyMeans.second - tmpPt.y) < mYDisTh){
                        curAreaCluster[k].emplace_back(tmpPt);
                        break;
                    }
                    // 如果是当前点云不符合聚类范围要求，且已经是最后一个已存在的聚类，那么以将该点建立一个新的聚类簇
                    else if(abs(distXY(tmpPt) - means) >= mXYDisTh && k == clusterCnt-1){
                        curAreaCluster.emplace_back(dummy);
                        curAreaCluster[clusterCnt][0] = tmpPt;
                    }
                    else{ 
                        continue; 
                    }
                    
                }
            }

            // 判断每个扇形区域中的聚类簇中点云数量是否符合要求
            int clusterCnt = curAreaCluster.size();
            for(int clusterNum = 0; clusterNum < clusterCnt; clusterNum++)
            {
                int ptCnt = curAreaCluster[clusterNum].size();
                if(ptCnt < 10){
                    continue;
                }
                tmpClustered.emplace_back(curAreaCluster[clusterNum]);
            }

        }

        //Merge the adjacent clusters 
        int clusterCnt = tmpClustered.size(); // 区域内所有的聚类簇数量
        
        std::vector<bool> toBeMerge(clusterCnt, false); // 标记该聚类簇是否已经合并
        std::multimap<int, int> mToBeMergeInd; // 记录可以合并的两个聚类簇ID
        std::set<int> sNeedMergeInd; // 需要合并的聚类簇中引导簇ID —— 因为后续是以该簇为标准，计算有哪些簇可以和该簇合并。用这个ID在mToBeMergeInd中搜索其余簇ID

        // Step 2 找到哪些簇需要合并
        for(int i = 0; i < clusterCnt; i++)
        {
            if(toBeMerge[i]){
                continue;
            }

            float means1 = computeClusterMeans(tmpClustered[i]);
            std::pair<float, float> xyMeans1;
            computeXYMeans(tmpClustered[i], xyMeans1);

            for(int j = 1; j < clusterCnt; j++)
            {
                if(toBeMerge[j]){
                    continue;
                }
                float means2 = computeClusterMeans(tmpClustered[j]);
                std::pair<float, float> xyMeans2;
                computeXYMeans(tmpClustered[j], xyMeans2);

                // 如果其它聚类簇中“均值点”距离很近，那么和当前
                if(abs(means1 - means2) < 0.6 && abs(xyMeans1.first - xyMeans2.first) < 0.8 && abs(xyMeans1.second - xyMeans2.second) < 0.8) 
                {
                    mToBeMergeInd.insert(std::make_pair(i, j));
                    sNeedMergeInd.insert(i);
                    toBeMerge[i] = true;
                    toBeMerge[j] = true;
                }
            }

        }

        // Step 3 合并聚类簇
        // 如果没有需要合并的簇，那么每个簇就是最后的聚类结果
        if(sNeedMergeInd.empty()){
            for(int i = 0; i < clusterCnt; i++)
            {
                clustered.emplace_back(tmpClustered[i]);
            }
        }
        else
        {
            // toBeMerge[i] == false说明该簇不需要和其它簇进行合并
            for(int i = 0; i < clusterCnt; i++)
            {
                if(toBeMerge[i] == false)
                {
                    clustered.emplace_back(tmpClustered[i]);
                }
            }
            
            // 
            for(auto setIt = sNeedMergeInd.begin(); setIt != sNeedMergeInd.end(); ++setIt)
            {
                // 需要和ID为needMergeInd的簇合并的簇数量有entries，且其ID记录在向量iter中
                int needMergeInd = *setIt;
                auto entries = mToBeMergeInd.count(needMergeInd);
                auto iter = mToBeMergeInd.find(needMergeInd);

                // 取出其余待合并簇的ID，记录在vInd
                std::vector<int> vInd; 
                while(entries){
                    int ind = iter->second;
                    vInd.emplace_back(ind);
                    ++iter;
                    --entries;
                }
                clustered.emplace_back(tmpClustered[needMergeInd]);
                size_t cnt = clustered.size();
                for(size_t j = 0; j < vInd.size(); j++)
                {
                    for(size_t ptNum = 0; ptNum < tmpClustered[vInd[j]].size(); ptNum++)
                    {
                        clustered[cnt - 1].emplace_back(tmpClustered[vInd[j]][ptNum]);
                    }
                }
            }
        }       

    }

    void LinK3D_Extractor::computeDirection(pcl::PointXYZ ptFrom, pcl::PointXYZ ptTo, Eigen::Vector2f &direction)
    {
        direction(0, 0) = ptTo.x - ptFrom.x;
        direction(1, 0) = ptTo.y - ptFrom.y;
    }

    /// @brief 将每个聚类簇的质心作为聚类关键点，并记录其对应的聚类簇ID
    /// @param clustered 
    /// @param index 返回关键点对应的聚类簇ID
    /// @return 
    std::vector<pcl::PointXYZ> LinK3D_Extractor::getAggregationKeyPt(const MatPt &clustered, std::vector<int> &index)
    {
        std::vector<pcl::PointXYZ> keyPoints;
        int clusterCnt = clustered.size();
        for(int clusterNum = 0; clusterNum < clusterCnt; clusterNum++)
        {
            int ptCnt = clustered[clusterNum].size();  
   
            if(ptCnt < mPtNumTh){
                continue;
            }

            std::vector<PointXYZSCA> tmpCluster;
            std::set<int> scans;
            float x = 0, y = 0, z = 0;
            for(int ptNum = 0; ptNum < ptCnt; ptNum++)
            {
                PointXYZSCA pt = clustered[clusterNum][ptNum];          
                int scan = pt.scan;
                scans.insert(scan);

                x += pt.x;
                y += pt.y;
                z += pt.z;
            }

            if((int)scans.size() < mScanNumTh){
                continue;
            }

            pcl::PointXYZ pt(x/ptCnt, y/ptCnt, z/ptCnt);
            keyPoints.emplace_back(pt);
            index.emplace_back(clusterNum); 
        }
        return keyPoints;
    }

    /// @brief 生成关键点描述子
    /// Step 1 计算所有keypoints中两点间的距离和方向
    /// Step 2 选择三个距离当前keypoint最近的关键点
    /// Step 3 对每一个keypoint都找到三个邻近点；分别以该关键点和其余邻近点的方向作为主方向，计算出三个不同的描述子。最后是180维的每个维度上选择距离最小的作为最终value
    /// @param keyPoints 
    /// @param descriptors 
    void LinK3D_Extractor::getDescriptor(const std::vector<pcl::PointXYZ> &keyPoints, cv::Mat &descriptors)
    {
        if(keyPoints.size() < 5){
            cout << "The number of keypoints is too small!" << endl;
            return;
        }

        size_t ptSize = keyPoints.size();

        // 一个描述子的维度是180维
        descriptors = cv::Mat::zeros(ptSize, 180, CV_32FC1); 

        // disTab和directionTab均是为了防止重复计算
        // disTab 对应论文中的 distance table 
        vector<vector<float>> disTab;
        vector<float> oneRowDis(ptSize, 0);
        disTab.resize(ptSize, oneRowDis);

        // directionTab 对应论文中的 direction table 
        vector<vector<Eigen::Vector2f>> directionTab;
        Eigen::Vector2f direct(0, 0);
        vector<Eigen::Vector2f> oneRowDirect(ptSize, direct);
        directionTab.resize(ptSize, oneRowDirect);

        // Step 1 计算所有keypoints中两点间的距离和方向
        for(size_t i = 0; i < keyPoints.size(); i++)
        {
            for(size_t j = i+1; j < keyPoints.size(); j++)
            {
                float dist = distPt2Pt(keyPoints[i], keyPoints[j]);
                disTab[i][j] = fRound(dist);
                disTab[j][i] = disTab[i][j];

                Eigen::Vector2f tmpDirection;
                                
                tmpDirection(0, 0) = keyPoints[j].x - keyPoints[i].x;
                tmpDirection(1, 0) = keyPoints[j].y - keyPoints[i].y;

                directionTab[i][j] = tmpDirection;
                directionTab[j][i] = -tmpDirection;
            }
        }

        
        for(size_t i = 0; i < keyPoints.size(); i++)
        {
            // Step 2 选择三个距离当前keypoint最近的关键点
            // disTab[i]中存储量其它距离keypoint i 的距离
            vector<float> tempRow(disTab[i]);
            std::sort(tempRow.begin(), tempRow.end()); // 距离升序排序

            //The number of the used closest keypoints is set to 3
            int Index[3];  
            for(int k = 0; k < 3; k++)  
            {
                vector<float>::iterator it1 = find(disTab[i].begin(), disTab[i].end(), tempRow[k+1]);
                if(it1 == disTab[i].end()){
                    continue;
                }else{
                    Index[k] = std::distance(disTab[i].begin(), it1);
                }
            }

            // Step 3 对每一个keypoint都找到三个邻近点；分别以该关键点和其余邻近点的方向作为主方向，计算出三个不同的描述子。最后是180维的每个维度上选择距离最小的作为最终value
            // 论文里也有详细的描述 —— 是为了降低描述子生成算法对邻近点的敏感度
            for(int indNum = 0; indNum < 3; indNum++) 
            {
                size_t index = Index[indNum];
                Eigen::Vector2f mainDirect;
                
                // 将当前keypoint和其距离最近点的方向作为主方向
                mainDirect = directionTab[i][index];

                //180-dimensional descriptor corresponds to 180 sector areas.
                vector<vector<float>> areaDis(180);  
                areaDis[0].emplace_back(disTab[i][index]);
                
                // 该for循环完成的工作可参见论文公式(3)
                // 将当前keypoint和其它点的方向作为XoY平面扇形区域划分的判断标准
                for(size_t j = 0; j < keyPoints.size(); j++)
                {
                    if(j == i || j == index){
                        continue;
                    }
                    
                    Eigen::Vector2f otherDirect = directionTab[i][j];

                    // 方向矩阵判断实际角度，因为acos()范围是[0, pi]。但是如果以Keypoint主方向为旋转轴，那么其它向量和主方向夹角应该在[0, 360]之间
                    // 用matrixDirect行列式正负判断：正值 [0, pi]， 负值[pi, 2*pi]
                    Eigen::Matrix2f matrixDirect;
                    matrixDirect << mainDirect(0, 0) , mainDirect(1, 0) , otherDirect(0, 0) , otherDirect(1, 0);
                    float deter = matrixDirect.determinant();

                    int areaNum = 0;
                    double cosAng = (double)mainDirect.dot(otherDirect) / (double)(mainDirect.norm() * otherDirect.norm());                                 
                    if(abs(cosAng) - 1 > 0){   
                        continue;
                    }
                    
                    float angel = acos(cosAng) * 180 / M_PI;
                    
                    if(angel < 0 || angel > 180){
                        continue;
                    }
                    
                    if(deter > 0){
                        areaNum = ceil((angel - 1) / 2); 
                        
                    }else{
                        if(angel - 2 < 0){ 
                            areaNum = 0;
                        }else{
                            angel = 360 - angel;
                            areaNum = ceil((angel - 1) / 2); 
                        }   
                    }

                    // !!!将两点之间的距离作为对应描述子维度上的value
                    if(areaNum != 0){
                        areaDis[areaNum].emplace_back(disTab[i][j]);
                    }
                            
                }
                
                //Each dimension of the descriptor corresponds to the distance between the current keypoint and the closest keypoint in the corresponding sector area.
                // 选择最小的距离值作为对应描述子维度上的value
                // TODO 这个循环放在整个for外面是否会好一点？
                float *descriptor = descriptors.ptr<float>(i);
                for(int areaNum = 0; areaNum < 180; areaNum++) 
                {
                    if(areaDis[areaNum].size() == 0){
                        continue;
                    }else{
                        std::sort(areaDis[areaNum].begin(), areaDis[areaNum].end());
                        if(descriptor[areaNum] == 0)
                        {
                            descriptor[areaNum] = areaDis[areaNum][0]; 
                        }                        
                    }
                }                
            }
        }        
    }

    float LinK3D_Extractor::fRound(float in)
    {
        float f;
        int temp = std::round(in * 10);
        f = temp/10.0;
        
        return f;
    }

    //The entrance of extracting LinK3D features
    /// @brief LinK3D特征点提取入口，最终会完成聚类关键点提取和对应描述子生成
    /// @param cloudIn 输入点云
    /// @param keyPoints 输出 聚类关键点
    /// @param descriptors 输出 聚类关键点对应的描述子
    /// @param index 输出 聚类关键点对应的聚类簇ID
    /// @param clustered 输出 聚类簇
    void LinK3D_Extractor::operator()(pcl::PointCloud<pcl::PointXYZ> &cloudIn, std::vector<pcl::PointXYZ> &keyPoints, cv::Mat &descriptors, std::vector<int> &index, MatPt &clustered)
    {
        if(cloudIn.empty()){
            return;
        }         

        MatPt scanCloud;
        // Step 1 提取边缘点
        getEdgePoint(cloudIn, scanCloud);
        
        MatPt areaCloud;
        // Step 2 扇形区域划分
        divideArea(scanCloud, areaCloud);       
        
        // Step 3 点云聚类
        computeCluster(areaCloud, clustered);
        
        // Step 4 获取聚类簇的质心，并将其作为聚类关键点以及对应的聚类簇ID
        keyPoints = getAggregationKeyPt(clustered, index);
        
        // Step 5 生成对应关键点的描述子
        getDescriptor(keyPoints, descriptors);         
    }

    /// @brief 根据聚类关键点完成两帧之间的匹配
    /// Step 1 两个for循环找到两帧中可能存在的匹配点对
    /// Step 2 在Step 1 中可能出现上一帧中的一个点对应当前帧中多个点，所以需要删除一对多的匹配
    /// @param descriptors1 输入 当前帧中聚类关键点描述子
    /// @param descriptors2 输入 上一帧中聚类关键点描述子
    /// @param vMatched 匹配关系
    void LinK3D_Extractor::matcher(cv::Mat &descriptors1, cv::Mat &descriptors2, vector<pair<int, int>> &vMatched)
    {
        int ptSize1 = descriptors1.rows;
        int ptSize2 = descriptors2.rows;
        
        std::multimap<int, int> matchedIndexScore;      
        std::multimap<int, int> mMatchedIndex;
        set<int> sIndex;
        
        // Step 1 两个for循环找到两帧中可能存在的匹配点对
        for(int i = 0; i < ptSize1; i++)
        {
            std::pair<int, int> highestIndexScore(0, 0); // 分别记录和当前descriptors1[1]匹配的descriptors2中的索引以及最高得分
            float* pDes1 = descriptors1.ptr<float>(i);        

            for(int j = 0; j < ptSize2; j++)
            {
                int sameBitScore = 0;
                float* pDes2 = descriptors2.ptr<float>(j);
                
                // descriptors1[i]和descriptors2[j]之间的得分计算
                for(int bitNum = 0; bitNum < 180; bitNum++)
                {
                    // 对应维度上的差值的绝对值 <= 0.2 那么相似度得分 +1
                    if(pDes1[bitNum] != 0 && pDes2[bitNum] != 0 && abs(pDes1[bitNum] - pDes2[bitNum]) <= 0.2){
                        sameBitScore += 1;
                    }             
                }

                if(sameBitScore > highestIndexScore.second){
                    highestIndexScore.first = j;
                    highestIndexScore.second = sameBitScore;
                }

            }
            //According to i, get the score
            matchedIndexScore.insert(std::make_pair(i, highestIndexScore.second)); 
            //According to j, get i
            mMatchedIndex.insert(std::make_pair(highestIndexScore.first, i)); 
            sIndex.insert(highestIndexScore.first);
        }
        
        //Remove one-to-multiple matches for descriptor2
        // Step 2 在Step 1 中可能出现上一帧中的一个点对应当前帧中多个点，所以需要删除一对多的匹配
        for(std::set<int>::iterator setIt = sIndex.begin(); setIt != sIndex.end(); ++setIt)
        {
            int indexJ = *setIt;
            auto entries = mMatchedIndex.count(indexJ);

            if(entries == 1){
                auto iterI = mMatchedIndex.find(indexJ);
                auto iterScore = matchedIndexScore.find(iterI->second);
                if(iterScore->second >= mScoreTh){   
                    vMatched.emplace_back(std::make_pair(iterI->second, indexJ));
                }           
            }
            // 出现一对多的情况，找到其中得分最高的点对作为最终匹配点对
            else{ 
                auto iter1 = mMatchedIndex.find(indexJ);
                int highestScore = 0;
                int highestScoreIndex = -1;

                while(entries){
                    int indexI = iter1->second;
                    auto iterScore = matchedIndexScore.find(indexI);
                    if(iterScore->second > highestScore){
                        highestScore = iterScore->second;
                        highestScoreIndex = indexI;
                    }                
                    ++iter1;
                    --entries;
                }
                if(highestScore >= mScoreTh){   
                    vMatched.emplace_back(std::make_pair(highestScoreIndex, indexJ));
                }            
            }
            
        }
    }

    //Get the true edge keypoints with higher smoothness
    /// @brief 获取具有高平滑度的边缘关键点
    /// 论文中有写到：不稳定的边缘点是分散的，有效的边缘点是垂直分布在簇中的，可见图 Fig. 4
    /*
        个人理解：并没有直接取确定Z轴的值后 统计当前位置上的关键点簇分布。
        具体做法是：通过各个聚类簇中的点在scan上的分布 —— 代替 Z轴位置，实现垂直方向上的划分
    */ 
    /// @param clustered 
    /// @param filtered 
    void LinK3D_Extractor::filterLowSmooth(MatPt &clustered, MatPt &filtered)
    {
        if(clustered.empty()){
            return;
        }

        int clusterSize = clustered.size();

        filtered.resize(clusterSize);
        // 遍历所有聚类簇
        for(int i = 0; i < clusterSize; i++)
        {
            int ptCnt = clustered[i].size();
            MatPt tmpCluster;
            vector<int> vScanID;
            // 遍历当前聚类簇中的所有点
            for(int j = 0; j < ptCnt; j++)
            {
                PointXYZSCA pt = clustered[i][j];
                int scan = pt.scan;
                auto it = std::find(vScanID.begin(), vScanID.end(), scan);
                if(it == vScanID.end()){
                    // vScanID.begin()的值实际上就是 j=0对应点的scan值，也就是当前聚类簇中第一个点的scan值
                    vScanID.emplace_back(scan);
                    vector<PointXYZSCA> vPt(1, pt);
                    tmpCluster.emplace_back(vPt);
                }else{
                    // filteredInd计算的是当前pt的scan值在vScanID中的位置距离第一个元素的距离
                    // 其实并没有考虑scan大小的排序。
                    /*
                        举个例子：比如vScanID.begin() = 3，后续遍历五个点的scan值分别是 4 1 8 3 2
                        那么这个时候vScanID的具体形式为 = [3, 4, 1, 8, 2]
                        后续遍历五个点的filteredInd分别为[1, 2, 3, 0, 4]
                        说明后续遍历的这五个点中的第四个点和j = 0的点被重新划分到一个簇中
                    */
                    // 所以，这个for的主要任务就是再根据scan值重新聚类每个点
                    int filteredInd = std::distance(vScanID.begin(), it);
                    tmpCluster[filteredInd].emplace_back(pt);
                }
            }

            // 从在同一条scan上的所有点中选出一个平滑度最高的点保留下来
            // ！！！如果一组垂直分布的边缘点分别位于不同scan上，可以大致认为该边缘点就是在这条scan上具有最高的平滑度。
            // 从这个思维反向出发 就能得到代码的目的
            for(size_t scanCnt = 0; scanCnt < tmpCluster.size(); scanCnt++)
            {
                if(tmpCluster[scanCnt].size() == 1){
                    filtered[i].emplace_back(tmpCluster[scanCnt][0]);
                }
                else
                {
                    float maxSmooth = 0;
                    PointXYZSCA maxSmoothPt;
                    for(size_t j = 0; j < tmpCluster[scanCnt].size(); j++)
                    {
                        if(tmpCluster[scanCnt][j].smoothness > maxSmooth){
                            maxSmooth = tmpCluster[scanCnt][j].smoothness;
                            maxSmoothPt = tmpCluster[scanCnt][j];
                        }
                    }
                    filtered[i].emplace_back(maxSmoothPt);
                }
            }  
        }
    }
    
    /// @brief 根据根据聚类关键点匹配结果，找到两帧间位于同一条scan值上的edge keypoint匹配点对，记录在matchedPoint
    /// @param filtered1 输入 当前帧上边缘关键点聚类簇
    /// @param filtered2 输入 上一帧上边缘关键点聚类簇
    /// @param index1 输入 当前帧上聚类簇的索引，便于根据索引找到对应的簇
    /// @param index2 
    /// @param matchedPoint 输出 帧间位于同一条scan值上的edge keypoint匹配点对
    /// @param vMatched 输入 聚类关键点匹配结果
    void LinK3D_Extractor::matchEdgePoints(MatPt &filtered1, MatPt &filtered2, vector<int> &index1, vector<int> &index2, MatPt &matchedPoint, vector<pair<int, int>> &vMatched)
    {
        if(vMatched.empty()){
            return;
        }

        for(size_t i = 0; i < vMatched.size(); i++)
        {
            // 从根据聚类关键点得到的匹配组中分别取出当前帧和上一帧中的ID关系
            pair<int, int> matchedInd = vMatched[i];
            int ind1 = index1[matchedInd.first];
            int ind2 = index2[matchedInd.second];
            
            int ptSize1 = filtered1[ind1].size();
            int ptSize2 = filtered2[ind2].size();

            // 以下两个map数据中 first = 簇中点对应的scan值；second = 簇中点的索引
            std::map<int, int> mScanID_Ind1;
            std::map<int, int> mScanID_Ind2;

            // mScanID_Ind1中记录 当前簇中所有点对应的scan和其索引
            for(int ptNum1 = 0; ptNum1 < ptSize1; ptNum1++)
            {
                int scanID1 = filtered1[ind1][ptNum1].scan;
                pair<int, int> scanID_Ind(scanID1, ptNum1);
                mScanID_Ind1.insert(scanID_Ind);
            }

            for(int ptNum2 = 0; ptNum2 < ptSize2; ptNum2++)
            {
                int scanID2 = filtered2[ind2][ptNum2].scan;
                pair<int, int> scanID_Ind(scanID2, ptNum2);
                mScanID_Ind2.insert(scanID_Ind);
            }

            // 找到两帧中位于同一scan上的egde keypoint的匹配对，记录在matchedPoint中
            for(auto it1 = mScanID_Ind1.begin(); it1 != mScanID_Ind1.end(); it1++)
            {
                int scanID1 = (*it1).first;
                auto it2 = mScanID_Ind2.find(scanID1);
                if(it2 == mScanID_Ind2.end()){
                    continue;
                }
                else
                {
                    vector<PointXYZSCA> tmpMatchPt;
                    PointXYZSCA pt1 = filtered1[ind1][(*it1).second];
                    PointXYZSCA pt2 = filtered2[ind2][(*it2).second];
                    tmpMatchPt.emplace_back(pt1);
                    tmpMatchPt.emplace_back(pt2);
                    matchedPoint.emplace_back(tmpMatchPt);
                }
            }
        }
    }

}

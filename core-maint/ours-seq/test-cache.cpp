#include <def.h>
void ClearCache()
{
    for(uint64 i = 0; i < tempRenderTargets.Count(); ++i)
    {
        TempRenderTarget* tempRT = tempRenderTargets[i];
        tempRT->RT.Shutdown();
        delete tempRT;
    }

    tempRenderTargets.RemoveAll(nullptr);

    for(uint64 i = 0; i < pipelineStates.Count(); ++i)
        DX12::DeferredRelease(pipelineStates[i].PSO);

    pipelineStates.RemoveAll(CachedPSO());
}

void main(){
    ClearCache();
}
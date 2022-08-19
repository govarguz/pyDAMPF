import asyncio
from motivus.client import Client

async def motivus_task(naux):
    motivus = await Client.connect()

    task_definition = {
        # local build
        'wasm_path': "build/mypmmoda-0.1.0.wasm",
        'loader_path': "build/mypmmoda-0.1.0.js",
        'data_link': "https://down-production.s3.amazonaws.com/mypmmoda-0.1.0.data.zip",
        # published build
        # 'algorithm': "mypmmoda",
        # 'algorithm_version': "0.1.0",
        'params': ['mypmmoda_MOD_mainbim', naux],
        'preload_files': {
            '/inputs.ini': motivus.read_file("./inputs.ini")
        },
        'result_files': {
            '/inputs.ini': './inputs.ini',
            '/zcdom.dfo': './zcdom.dfo',
            '/tdom.dfo': './tdom.dfo'
        }
    }

    task_id = motivus.call_async(task_definition)
    task = motivus.select_task(task_id)
    return await task

class mypmmoda(object):
    @classmethod
    def mainbim(cls, naux):
        asyncio.run(motivus_task(naux))

if __name__ == '__main__':
    mypmmoda.mainbim(1)

